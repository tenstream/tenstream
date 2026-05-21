#include <stdalign.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

/* Storage is provided by Fortran bind(C) module variables in c_syscall_wrappers.F90.
 * This constructor initialises them from system headers before any Fortran code runs. */
extern int c_SC_PAGESIZE;
extern int cF_ULOCK;
extern int cF_LOCK;
extern int cPROT_READ;
extern int cMAP_PRIVATE;
extern int cMAP_NORESERVE;
extern int cO_RDONLY;
extern int cO_CREAT;
extern int cO_RDWR;
extern int cO_WRONLY;
extern int cO_APPEND;
extern int cS_IRWXU;
extern int cS_IRUSR;
extern int cS_IWUSR;
extern int cS_IROTH;
extern int default_user_wrmode;

__attribute__((constructor))
static void init_c_syscall_constants(void) {
  c_SC_PAGESIZE    = (int)_SC_PAGESIZE;
  cF_ULOCK         = (int)F_ULOCK;
  cF_LOCK          = (int)F_LOCK;
  cPROT_READ       = (int)PROT_READ;
  cMAP_PRIVATE     = (int)MAP_PRIVATE;
  cMAP_NORESERVE   = (int)MAP_NORESERVE;
  cO_RDONLY        = (int)O_RDONLY;
  cO_CREAT         = (int)O_CREAT;
  cO_RDWR          = (int)O_RDWR;
  cO_WRONLY        = (int)O_WRONLY;
  cO_APPEND        = (int)O_APPEND;
  cS_IRWXU         = (int)S_IRWXU;
  cS_IRUSR         = (int)S_IRUSR;
  cS_IWUSR         = (int)S_IWUSR;
  cS_IROTH         = (int)S_IROTH;
  default_user_wrmode = (int)(S_IRUSR | S_IWUSR);
}
