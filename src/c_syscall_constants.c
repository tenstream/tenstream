#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

const int cF_ULOCK = F_ULOCK;
const int cF_LOCK = F_LOCK;

const int c_SC_PAGESIZE = _SC_PAGESIZE;

const int cPROT_READ      = PROT_READ;
const int cMAP_PRIVATE    = MAP_PRIVATE;
const int cMAP_NORESERVE  = MAP_NORESERVE;
const int cO_RDONLY       = O_RDONLY;
const int cO_CREAT        = O_CREAT;
const int cO_RDWR         = O_RDWR;
const int cO_WRONLY       = O_WRONLY;
const int cO_APPEND       = O_APPEND;
const int cS_IRWXU        = S_IRWXU;
const int cS_IRUSR        = S_IRUSR;
const int cS_IWUSR        = S_IWUSR;
const int cS_IROTH        = S_IROTH;

const mode_t default_user_wrmode = S_IRUSR|S_IWUSR;
