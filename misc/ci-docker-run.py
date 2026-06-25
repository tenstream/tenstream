#!/usr/bin/env python3
"""
ci-docker-run.py  —  Run a .gitlab-ci.yml job locally in Docker.

Usage:
    ./misc/ci-docker-run.py                  # interactive menu
    ./misc/ci-docker-run.py <job-name>       # run a specific job
    ./misc/ci-docker-run.py --list           # list available jobs
    ./misc/ci-docker-run.py --dry-run <job>  # print the docker command without running

Options:
    --all           Include all stages (default: only build_docker stage)
    --no-tty        Disable TTY allocation (-t flag) for docker run
    --dry-run       Print the generated script and docker command, don't run
    --list          List available jobs and exit
"""

import os
import sys
import argparse
import subprocess
import tempfile


# ── YAML loading ──────────────────────────────────────────────────────────────

def load_yaml(path):
    try:
        import yaml
    except ImportError:
        print("PyYAML not found — installing...", flush=True)
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--quiet", "pyyaml"]
        )
        import yaml  # noqa: PLC0415
    with open(path) as fh:
        return yaml.safe_load(fh)


# ── Job resolution ────────────────────────────────────────────────────────────

def _merge(base, override):
    """Deep-merge two dicts. Dicts recurse; everything else: override wins."""
    if not isinstance(base, dict) or not isinstance(override, dict):
        return override
    result = dict(base)
    for k, v in override.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _merge(result[k], v)
        else:
            result[k] = v
    return result


def resolve(name, jobs, _seen=None):
    """Recursively resolve a job's config by following extends chains."""
    if _seen is None:
        _seen = set()
    if name in _seen:
        raise ValueError(f"Circular extends detected: {name}")
    _seen = _seen | {name}

    raw = dict(jobs[name])
    extends = raw.pop("extends", None)
    if extends is None:
        return raw

    if isinstance(extends, str):
        extends = [extends]

    merged_parent = {}
    for parent_name in extends:
        if parent_name not in jobs:
            raise ValueError(f"Job '{name}' extends unknown job '{parent_name}'")
        merged_parent = _merge(merged_parent, resolve(parent_name, jobs, _seen))

    return _merge(merged_parent, raw)


# ── Helpers ───────────────────────────────────────────────────────────────────

_NON_JOB_KEYS = {"stages", "variables", "image", "before_script", "after_script",
                 "include", "workflow", "default", "cache"}


def flatten(lst):
    """Flatten one level of nesting — needed because PyYAML expands YAML anchors
    that reference a sequence into a nested list when used as a sequence item."""
    out = []
    for item in lst:
        if isinstance(item, list):
            out.extend(item)
        else:
            out.append(item)
    return out


def is_runnable(name, job):
    return (
        not name.startswith(".")
        and isinstance(job, dict)
        and name not in _NON_JOB_KEYS
        and "script" in job
    )


def image_name(image):
    if isinstance(image, dict):
        return image.get("name", "ubuntu:latest")
    return image or "ubuntu:latest"


def as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return flatten(value)
    return [value]


def expand_vars(variables):
    """One-pass expansion of ${VAR}/$VAR references within the variable dict,
    mimicking GitLab CI's own variable interpolation at pipeline time."""
    import re
    result = dict(variables)
    for key in list(result):
        val = str(result[key])
        def repl(m):
            name = m.group(1) or m.group(2)
            return str(result.get(name, m.group(0)))
        result[key] = re.sub(r'\$\{(\w+)\}|\$(\w+)', repl, val)
    return result


# ── Core ──────────────────────────────────────────────────────────────────────

def build_inner_script(job, global_before, ci_project_dir):
    lines = ["set -e", "set -o pipefail", "set -x",
             "trap 'echo; echo \"--- command failed (exit $?) --- dropping to shell\"; exec bash -i' ERR",
             ""]

    before = global_before + as_list(job.get("before_script"))
    script = as_list(job.get("script"))

    if before:
        lines.append("# --- before_script ---")
        lines.extend(before)
        lines.append("")

    lines.append("# --- script ---")
    lines.extend(script)
    return "\n".join(lines)


def run_job(job_name, all_jobs, global_vars, global_before,
            repo_root, dry_run=False, tty=True):
    job = resolve(job_name, all_jobs)

    img = image_name(job.get("image"))

    variables = dict(global_vars)
    variables.update(job.get("variables") or {})

    ci_project_dir = "/builds/tenstream/tenstream"
    variables.setdefault("CI_PROJECT_DIR", ci_project_dir)
    variables.setdefault("DEBIAN_FRONTEND", "noninteractive")

    variables = expand_vars(variables)

    # Git 2.35.2+ refuses to operate on repos owned by a different user.
    # Inside Docker the process runs as root but mounted files belong to the
    # host user, triggering the "dubious ownership" check. Pass the exception
    # via env vars so it applies even before git is installed.
    variables["GIT_CONFIG_COUNT"] = "1"
    variables["GIT_CONFIG_KEY_0"] = "safe.directory"
    variables["GIT_CONFIG_VALUE_0"] = "*"

    inner = build_inner_script(job, global_before, ci_project_dir)

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".sh", delete=False, prefix="ci_job_"
    )
    tmp.write(inner)
    tmp.close()

    container_script = "/tmp/ci_job_run.sh"

    docker_cmd = ["docker", "run", "--rm"]
    if tty:
        docker_cmd.append("-it")
    else:
        docker_cmd.append("-i")

    docker_cmd += ["--workdir", ci_project_dir]
    docker_cmd += ["-v", f"{repo_root}:{ci_project_dir}"]
    docker_cmd += ["-v", f"{tmp.name}:{container_script}:ro"]

    for k, v in variables.items():
        docker_cmd += ["-e", f"{k}={v}"]

    docker_cmd += [img, "bash", container_script]

    print(f"\n{'='*64}")
    print(f"  Job   : {job_name}")
    print(f"  Image : {img}")
    if variables:
        print("  Env   :")
        for k, v in sorted(variables.items()):
            print(f"    {k}={v}")
    print(f"{'='*64}")

    if dry_run:
        print("\n--- inner script ---")
        for i, line in enumerate(inner.splitlines(), 1):
            print(f"  {i:3d}  {line}")
        print("\n--- docker command ---")
        print("  " + " \\\n    ".join(docker_cmd))
        print()
        os.unlink(tmp.name)
        return

    print(f"\nStarting container...\n")
    try:
        ret = subprocess.run(docker_cmd)
        sys.exit(ret.returncode)
    finally:
        os.unlink(tmp.name)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run a .gitlab-ci.yml job locally in Docker.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("job", nargs="?", help="Job name to run")
    parser.add_argument("--list", action="store_true", help="List available jobs and exit")
    parser.add_argument("--all", action="store_true",
                        help="Include all stages, not just build_docker")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print generated script and docker command without running")
    parser.add_argument("--no-tty", action="store_true",
                        help="Disable TTY allocation for docker run")
    args = parser.parse_args()

    # CI YAML lives one directory up from misc/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    ci_yaml = os.path.join(repo_root, ".gitlab-ci.yml")

    if not os.path.exists(ci_yaml):
        sys.exit(f"Error: {ci_yaml} not found")

    data = load_yaml(ci_yaml)

    all_jobs = {
        k: v for k, v in data.items()
        if isinstance(v, dict) and k not in _NON_JOB_KEYS
    }

    global_vars = dict(data.get("variables") or {})
    global_before = as_list(data.get("before_script"))

    concrete = {
        name: job for name, job in all_jobs.items()
        if is_runnable(name, resolve(name, all_jobs))
    }

    def is_docker_stage(name):
        job = resolve(name, all_jobs)
        return job.get("stage", "") == "build_docker"

    display = concrete if args.all else {
        k: v for k, v in concrete.items() if is_docker_stage(k)
    }

    if args.list:
        print("Available jobs:")
        for name in sorted(display):
            job = resolve(name, all_jobs)
            img = image_name(job.get("image", ""))
            stage = job.get("stage", "")
            print(f"  {name:<40s}  [{img}]  (stage: {stage})")
        return

    job_name = args.job

    if job_name is None:
        names = sorted(display)
        if not names:
            sys.exit("No runnable jobs found in .gitlab-ci.yml")
        print("Available CI jobs:")
        for i, n in enumerate(names, 1):
            job = resolve(n, all_jobs)
            img = image_name(job.get("image", ""))
            print(f"  {i:2d}.  {n:<40s}  [{img}]")
        print()
        choice = input("Select job number or name: ").strip()
        if choice.isdigit():
            idx = int(choice) - 1
            if not (0 <= idx < len(names)):
                sys.exit("Invalid selection")
            job_name = names[idx]
        else:
            job_name = choice

    if job_name not in all_jobs:
        sys.exit(
            f"Unknown job '{job_name}'.\n"
            f"Run with --list to see available jobs."
        )

    run_job(
        job_name,
        all_jobs,
        global_vars,
        global_before,
        repo_root=repo_root,
        dry_run=args.dry_run,
        tty=not args.no_tty,
    )


if __name__ == "__main__":
    main()
