#include <stdalign.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

alignas(8) const int cF_ULOCK = F_ULOCK;
alignas(8) const int cF_LOCK = F_LOCK;

alignas(8) const int c_SC_PAGESIZE = _SC_PAGESIZE;

alignas(8) const int cPROT_READ      = PROT_READ;
alignas(8) const int cMAP_PRIVATE    = MAP_PRIVATE;
alignas(8) const int cMAP_NORESERVE  = MAP_NORESERVE;
alignas(8) const int cO_RDONLY       = O_RDONLY;
alignas(8) const int cO_CREAT        = O_CREAT;
alignas(8) const int cO_RDWR         = O_RDWR;
alignas(8) const int cO_WRONLY       = O_WRONLY;
alignas(8) const int cO_APPEND       = O_APPEND;
alignas(8) const int cS_IRWXU        = S_IRWXU;
alignas(8) const int cS_IRUSR        = S_IRUSR;
alignas(8) const int cS_IWUSR        = S_IWUSR;
alignas(8) const int cS_IROTH        = S_IROTH;

alignas(8) const mode_t default_user_wrmode = S_IRUSR|S_IWUSR;
