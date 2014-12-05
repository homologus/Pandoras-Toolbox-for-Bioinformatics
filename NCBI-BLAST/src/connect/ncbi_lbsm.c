/* $Id: ncbi_lbsm.c 440972 2014-07-18 16:07:40Z lavr $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   LBSM client-server data exchange API
 *
 *   UNIX only !!!
 *
 */

#include "ncbi_ansi_ext.h"
#include "ncbi_lbsm.h"
#include "ncbi_priv.h"
#include <connect/ncbi_socket_unix.h>
#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef   max
#  undef max
#endif /*max*/
#define  max(a, b)  ((a) < (b) ? (b) : (a))

#ifdef   min
#  undef min
#endif /*min*/
#define  min(a, b)  ((a) > (b) ? (b) : (a))

#ifdef   fabs
#  undef fabs
#endif /*fabs*/
#define  fabs(v)    ((v) < 0.0 ? -(v) : (v))

#define NCBI_USE_ERRCODE_X   Connect_LBSM

/* Alias for service lookups by type */
#define eLBSM_PendingOrService eLBSM_Invalid


int/*bool*/ LBSM_SetVersion(HEAP heap, const SLBSM_Version* v)
{
    SLBSM_Version* b;
    size_t size;

    if (!v  ||  v->entry.type != eLBSM_Version)
        return 0/*false: bad version*/;

    if (HEAP_Next(heap, 0))
        return 0/*false: heap is not empty*/;

    size = sizeof(*v) - sizeof(v->entry.head);
    if (!(b = (SLBSM_Version*) HEAP_Alloc(heap, size, 0))
        ||  (void*) b != (void*) HEAP_Base(heap)) {
        return 0/*false: bad allocation*/;
    }

    /* Copy version to the heap, preserving block header */
    memcpy((char*) b + sizeof(b->entry.head),
           (char*) v + sizeof(v->entry.head), size);
    return 1/*true*/;
}


const SLBSM_Version* LBSM_GetVersion(HEAP heap)
{
    const SLBSM_Entry* e = (const SLBSM_Entry*) HEAP_Next(heap, 0);
    if (!e  ||  e->type != eLBSM_Version)
        return 0;
    assert((void*) e == (void*) HEAP_Base(heap));
    return (const SLBSM_Version*) e;
}


int/*bool*/ LBSM_PutConfig(HEAP heap, const char* name)
{
    SLBSM_Entry* e = 0;
    SLBSM_Config* c;
    size_t size;

    while ((e = (SLBSM_Entry*) HEAP_Next(heap, &e->head)) != 0) {
        if (e->type == eLBSM_Config)
            HEAP_Free(heap, &e->head);
    }
    if (!name)
        name = "";

    size = sizeof(*c) - sizeof(c->entry.head);
    if (!(c = (SLBSM_Config*) HEAP_Alloc(heap, size + strlen(name), 0)))
        return 0/*failure*/;

    c->entry.good = NCBI_TIME_INFINITE;
    c->entry.type = eLBSM_Config;
    strcpy(c->name, name);

    return 1/*success*/;
}


const char* LBSM_GetConfig(HEAP heap)
{
    const SLBSM_Entry* e = 0;
    while ((e = (const SLBSM_Entry*) HEAP_Next(heap, &e->head)) != 0) {
        if (e->type == eLBSM_Config) {
            const SLBSM_Config* c = (const SLBSM_Config*) e;
            return c->name;
        }
    }
    return 0;
}


#ifdef __GNUC__
inline
#endif /*__GNUC__*/
static int/*bool*/ x_NameMatch(int/*bool*/ mask,
                               const char* name,
                               const char* test)
{
    return ((!mask  &&  strcasecmp      (name, test) == 0)  ||
            ( mask  &&  UTIL_MatchesMask(name, test)) ? 1/*T*/ : 0/*F*/);
}


static SLBSM_Service* s_Lookup(HEAP heap, const char* name, int/*bool*/ mask,
                               ELBSM_Type type, const SLBSM_Service* s)
{
    SLBSM_Entry* e = (SLBSM_Entry*) &s->entry;

    while ((e = (SLBSM_Entry*) HEAP_Next(heap, &e->head)) != 0) {
        if (e->type != eLBSM_Service  &&  e->type != eLBSM_Pending)
            continue/*not a server entry*/;
        if (type == eLBSM_PendingOrService  ||  type == e->type) {
            SLBSM_Service* t = (SLBSM_Service*) e;
            assert(t->info.host);
            if (!name  ||  x_NameMatch(mask, (const char*) t + t->name, name))
                return t;
        }
    }
    return 0;
}


const SLBSM_Service* LBSM_LookupService(HEAP                 heap,
                                        const char*          name,
                                        int/*bool*/          mask,
                                        const SLBSM_Service* prev)
{
    if (prev  &&  prev->entry.type != eLBSM_Service) {
        CORE_LOG_X(1, eLOG_Error, "Invalid preceding entry in service lookup");
        return 0;
    }
    return s_Lookup(heap, name, mask, eLBSM_Service, prev);
}


/* NB: "hint" may point to a freed and merged entry */
const SLBSM_Host* LBSM_LookupHost(HEAP               heap,
                                  unsigned int       addr,
                                  const SLBSM_Entry* hint)
{
    int/*bool*/ wrap = hint ? 1/*true*/ : 0/*false*/;
    const SLBSM_Entry* e = hint;
    while ((e = (const SLBSM_Entry*) HEAP_Next(heap, &e->head)) != hint) {
        if (!e) {
            assert(hint);
            if (!wrap)
                break;
            wrap = 0/*false*/;
            continue;
        }
        if (e->type == eLBSM_Host) {
            const SLBSM_Host* h = (const SLBSM_Host*) e;
            assert(h->addr);
            if (!addr  ||  h->addr == addr)
                return h;
        }
    }
    return 0;
}


ELBSM_Type LBSM_PublishService(HEAP heap, const SLBSM_Service* s, int log)
{
    SLBSM_Service* t = 0;
    SLBSM_Service* b;
    const char* name;
    ELBSM_Type type;
    size_t size;

    if (!s  ||  s->entry.type != eLBSM_Service  ||  !s->info.host)
        return eLBSM_Invalid/*bad call*/;

    name = (const char*) s + s->name;
    type = s->info.quorum ? eLBSM_PendingOrService : eLBSM_Service;
    while ((t = s_Lookup(heap, name, 0, type, t)) != 0) {
        char addr1[40], addr2[40];
        int/*bool*/ equal;

        if (s->info.type != t->info.type  ||
            s->info.host != t->info.host  ||
            s->info.port != t->info.port) {
            continue/*completely distinct*/;
        }
        if (!s->info.quorum != !t->info.quorum)
            continue/*backup does not preempt service and vice versa*/;
        if (s->info.quorum /*implies "t->info.quorum"*/) {
            type = t->entry.type;  /* can be both pending or already serving */
            break/*replace backup entry*/;
        }
        assert(type == eLBSM_Service);
        if (s->info.host == s->addr /*NB: implies "t->info.host == s->addr"*/)
            break/*replacement by the owner of the service*/;
        equal = SERV_EqualInfo(&s->info, &t->info);
        if (equal
            &&  (s->addr == t->addr  ||  s->info.rate * t->info.rate >= 0.0)
            &&  s->info.sful == t->info.sful
            &&  s->info.locl == t->info.locl
            &&  s->info.flag == t->info.flag) {
            if (s->info.host != t->addr/*unowned service*/
                &&  (!s->info.rate ^ !t->info.rate  ||
                     !s->info.time ^ !t->info.time)) {
                /* Unowned service is changing its state, so to avoid excessive
                 * flapping make sure state is sticky for at least half-life */
                if (s->entry.good <
                    t->entry.good + (max(t->info.time, s->info.time) >> 1)) {
                    return type/*fake success*/;
                }
            }
            break/*even definitions (by possibly 2 different hosts)*/;
        }
        if (log  &&  SOCK_ntoa(s->addr, addr1, sizeof(addr1)) != 0)
            sprintf(addr1, "%08X", s->addr);
        if (log  &&  SOCK_ntoa(t->addr, addr2, sizeof(addr2)) != 0)
            sprintf(addr2, "%08X", t->addr);
        if (equal) {
            /*even definitions but with different dynamic properties*/
            if (log) {
                CORE_LOGF_X(2, eLOG_Warning,
                            ("Service `%s' defined by"
                             " both %s and %s", name, addr1, addr2));
            }
            break;
        }
        /* Here we have UNEVEN definitions, and the new one is not run
         * by the announcing host, i.e. s->info.host != s->addr */
        if (s->info.host == t->addr /*NB: implies "t->info.host == t->addr"*/){
            if (log) {
                CORE_LOGF_X(3, eLOG_Warning,
                            ("Rejected attempt from %s to unevenly replace"
                             " service `%s' from %s", addr1, name, addr2));
            }
            return type/*fake success*/;
        }
        /* Here both hosts supply uneven definitions of the unowned service */
        if (log) {
            CORE_LOGF_X(4, eLOG_Error,
                        ("Service `%s' announced unevenly by"
                         " both %s and %s", name, addr1, addr2));
        }
        break;
    }
    if (t)
        HEAP_Free(heap, &t->entry.head);
    else if (!type)
        type = eLBSM_Pending;

    size  = (size_t)((name + strlen(name) + 1) - (char*) s);
    size -= sizeof(s->entry.head);
    if (!(b = (SLBSM_Service*) HEAP_Alloc(heap, size, 0)))
        return eLBSM_Invalid/*failure*/;

    /* Copy service to the heap, preserving block header */
    memcpy((char*) b + sizeof(b->entry.head),
           (char*) s + sizeof(s->entry.head), size);
    b->entry.type = type;
    return type/*success*/;
}


TNCBI_Time LBSM_UnpublishHost(HEAP heap, unsigned int addr, unsigned int items)
{
    SLBSM_Entry* p, *e = 0;
    TNCBI_Time joined = 0;

    /* Delete a host and [almost] all service entries announced by that host */
    for (p = 0;  (e = (SLBSM_Entry*) HEAP_Walk(heap, &e->head)) != 0;  p = e) {
        if (!(e->head.flag & 1))
            continue/*free block*/;
        if (e->type == eLBSM_Host) {
            SLBSM_Host* h = (SLBSM_Host*) e;
            assert(h->addr);
            if (h->addr != addr)
                continue/*not this host*/;
            assert(!joined);
            joined
                = h->sys.data.joined ? h->sys.data.joined : NCBI_TIME_INFINITE;
            --items;
        } else if (e->type == eLBSM_Service  ||  e->type == eLBSM_Pending) {
            SLBSM_Service* s = (SLBSM_Service*) e;
            assert(s->info.host);
            if (s->addr != addr)
                continue/*not by this host*/;
            --items;
            if (s->entry.type == eLBSM_Service
                &&  (s->info.rate < 0.0  ||  s->info.quorum)) {
                continue/*keep static / active backup entry*/;
            }
        } else
            continue;
        HEAP_FreeFast(heap, &e->head, &p->head);
        if (!items)
            break;
        if (p  &&  !(p->head.flag & 1))
            e = p;
    }
    return joined;
}


int/*bool*/ LBSM_PublishHost(HEAP heap, const SLBSM_Host* h)
{
    SLBSM_Host* b;
    size_t size;

    if (!h  ||  h->entry.type != eLBSM_Host  ||  !h->addr)
        return 0/*failure*/;

#if defined(_DEBUG)  &&  !defined(NDEBUG)
    /* This used to be an assert(!LBSM_LookupHost(heap, h->addr, 0)); */
    b = (SLBSM_Host*) LBSM_LookupHost(heap, h->addr, 0);
    if (b) {
        char addr[40];
        if (SOCK_ntoa(h->addr, addr, sizeof(addr)) != 0)
            sprintf(addr, "%08X", h->addr);
        CORE_LOGF(eLOG_Fatal, ("Host %s found at %p", addr, b));
    }
#endif /*_DEBUG && !NDEBUG*/

    /* Then, allocate a new host entry */
    size  = sizeof(*h) + (h->env ? strlen((const char*) h + h->env) + 1 : 0);
    size -= sizeof(h->entry.head);
    if (!(b = (SLBSM_Host*) HEAP_Alloc(heap, size, 1)))
        return 0/*failure*/;

    /* Copy host to the heap, preserving the block header */
    memcpy((char*) b + sizeof(b->entry.head),
           (char*) h + sizeof(h->entry.head), size);

    return 1/*success*/;
}


size_t LBSM_Expire(HEAP heap, TNCBI_Time time, size_t count)
{
    int/*bool*/ need_backup_rescan = 0;
    SLBSM_Entry* p, *e = 0;
    size_t n = 0;

    for (p = 0;  (e = (SLBSM_Entry*) HEAP_Walk(heap, &e->head)) != 0;  p = e) {
        if (!(e->head.flag & 1)           /* vacant block -- nothing to do */
            ||  e->type == eLBSM_Version  /* Version entry never expires   */
            ||  e->type == eLBSM_Config   /* Config entry never expires    */
            ||  e->good >= time) {        /* Entry is not yet expired      */
            continue;
        }
        if (e->type == eLBSM_Service) {
            const SLBSM_Service* s = (const SLBSM_Service*) e;
            assert(s->info.host);
            if (s->info.quorum) {
                e->good = 0;
                need_backup_rescan = 1;
                continue/*do not delete active backup entries just
                          as yet, but later at rescan time*/;
            }
        } else if (e->type == eLBSM_Host) {
            assert(((const SLBSM_Host*) e)->addr);
            if (count) {
                const SLBSM_Host* h = (const SLBSM_Host*) e;
                char addr[40];
                char buf[32];
                if (SOCK_ntoa(h->addr, addr, sizeof(addr)) != 0)
                    sprintf(addr, "%08X", h->addr);
                if (count != (size_t)(-1L)) {
                    sprintf(buf, " (%lu)", (unsigned long) --count);
                    if (!count)
                        count = (size_t)(-1L);
                } else
                    *buf = '\0';
                CORE_LOGF_X(5, eLOG_Warning, ("Host %s expired%s", addr, buf));
            }
            n++;
        }
        HEAP_FreeFast(heap, &e->head, &p->head);
        if (p  &&  !(p->head.flag & 1))
            e = p;
    }
    if (need_backup_rescan)
        LBSM_BackupWatch(heap, count ? 1 : 0);
    return n;
}


/*FIXME: Must be hard optimized!*/
static int s_BackupWatchByName(HEAP heap, const char* name, int/*bool*/ log)
{
    int/*bool*/ backed_up = 0;
    double host_status = 0.0;
    double serv_status = 0.0;
    SLBSM_Service* s = 0;
    int n_pending = 0;
    int n_active = 0;
    int quorum = 0;

    assert(sizeof(s->info.reserved) >= sizeof(double));
    while ((s = s_Lookup(heap, name, 0, eLBSM_PendingOrService, s)) != 0) {
        if (s->info.quorum) {
            const SLBSM_Service* t = 0;
            if (s->entry.type == eLBSM_Service) {
                s->entry.type  = eLBSM_Pending;
                backed_up = 1;
            }
            if (!s->entry.good) {
                HEAP_Free(heap, &s->entry.head);
                continue/*expired entry left behind by LBSM_Expire()*/;
            }
            if (!s->info.rate)
                continue;
            if (quorum > s->info.quorum  ||  !quorum)
                quorum = s->info.quorum;
            while ((t = s_Lookup(heap, name, 0, eLBSM_Service, t)) != 0) {
                if (t->info.type == s->info.type  &&
                    t->info.host == s->info.host  &&
                    t->info.port == s->info.port) {
                    break/*exact match in real services*/;
                }
            }
            s->info.reserved[sizeof(s->info.reserved) - 1] = !t ? 1 : 0;
            if (!t)
                n_pending++; /*unique*/
        } else if (s->info.rate) {
            double rate;
            memcpy(&rate, s->info.reserved + sizeof(s->info.reserved)
                   - sizeof(rate), sizeof(rate));
            if (rate) {
                assert(rate > 0.0);
                memset(s->info.reserved + sizeof(s->info.reserved)
                       - sizeof(rate), 0, sizeof(rate));
                s->info.rate = rate;
            }
            if (s->info.rate > 0.0) {
                const SLBSM_Host* h = LBSM_LookupHost(heap, s->info.host, 0);
                if (h) {
                    host_status += (s->info.flag & eSERV_Blast
                                    ? h->sys.load.statusBLAST
                                    : h->sys.load.status);
                } else
                    host_status += 1.0 / LBSM_DEFAULT_RATE;
                serv_status += s->info.rate;
            }
            n_active++;
        }
    }

    if ((!backed_up  &&   n_active <  quorum  &&   n_pending)  ||
        ( backed_up  &&  (n_active >= quorum  ||  !n_pending))) {
        backed_up = !backed_up;
        if (log) {
            CORE_LOGF_X(6, backed_up || n_active < quorum ? eLOG_Warning
                                                          : eLOG_Note,
                        ("%s `%s' [found %d, quorum %d]", backed_up
                         ? "Backup for" : n_active < quorum
                         ? "No backup for " : "Recovered",
                         name, n_active, quorum));
        }
    }
    if (!backed_up)
        return 0/*false*/;

    assert(s == 0);
    while ((s = s_Lookup(heap, name, 0, eLBSM_PendingOrService, s)) != 0) {
        if (!s->info.rate)
            continue;
        if (s->entry.type == eLBSM_Service) {
            if (s->info.quorum)
                continue;
            if (s->info.rate > 0.0) {
                const SLBSM_Host* h = LBSM_LookupHost(heap, s->info.host, 0);
                double rate = s->info.rate;
                memcpy(s->info.reserved + sizeof(s->info.reserved)
                       - sizeof(rate), &rate, sizeof(rate));
                if (h) {
                    rate = (s->info.flag & eSERV_Blast
                            ? h->sys.load.statusBLAST
                            : h->sys.load.status);
                } else
                    rate = 1.0 / LBSM_DEFAULT_RATE;
                s->info.rate = -rate * serv_status / host_status;
            }
        } else {
            assert(s->entry.type == eLBSM_Pending  &&  s->info.quorum);
            if (!s->info.reserved[sizeof(s->info.reserved) - 1]
                ||  n_active >= quorum) {
                continue;
            }
            do {
                SLBSM_Service* q = s;
                SLBSM_Service* t = q;
                while ((t = s_Lookup(heap, name, 0, eLBSM_Pending, t)) != 0) {
                    if (t->info.reserved[sizeof(t->info.reserved) - 1]
                        &&  t->info.quorum < q->info.quorum)
                        q = t;
                }
                q->info.reserved[sizeof(q->info.reserved) - 1] = 0;
                q->entry.type = eLBSM_Service;
                n_active++;
                if (q == s)
                    break;
            } while (n_active < quorum);
        }
    }
    return 1/*true*/;
}


int/*bool*/ LBSM_BackupWatch(HEAP heap, int/*bool*/ log)
{
    int backed_up = 0/*false*/;
    const SLBSM_Service* s = 0;

    while ((s = s_Lookup(heap, 0/*any*/, 0, eLBSM_PendingOrService, s)) != 0) {
        if (s->info.quorum) {
            const char* name = (const char*) s + s->name;
            const SLBSM_Service* t = 0;
            while ((t = s_Lookup(heap,name,0,eLBSM_PendingOrService,t)) != 0) {
                if (t->info.quorum)
                    break;
            }
            if (s == t  &&  s_BackupWatchByName(heap, name, log))
                backed_up = 1/*true*/;
        }
    }
    return backed_up;
}


double LBSM_CalculateStatus(double rate, double fine, ESERV_Flag flag,
                            const SLBSM_HostLoad* load)
{
    double status;

    if (!rate)
        return 0.0;
    if (rate < LBSM_STANDBY_THRESHOLD)
        status = rate < 0.0 ? -LBSM_DEFAULT_RATE : LBSM_DEFAULT_RATE;
    else
        status = flag & eSERV_Blast ? load->statusBLAST : load->status;
    status *= rate / LBSM_DEFAULT_RATE;
    /* accurately apply fine: avoid any imperfections with 100% */
    status *= (100. - (fine < 0. ? 0. : fine > 100. ? 100. : fine)) / 100.0;
    return fabs(status); /*paranoid but safe*/
}


/* NB: This call specifically omits access control (uid) so that the API
 * can be used from the dispatcher (user www) [cf. command line tool].
 */
int/*bool*/ LBSM_SubmitPenaltyOrRerate(const char*    name,
                                       ESERV_Type     type,
                                       double         rate,
                                       TNCBI_Time     fine,
                                       unsigned int   host,
                                       unsigned short port,
                                       const char*    path)
{
    const char* type_name = type ? SERV_TypeStr(type) : "ANY";
    struct sigaction sa, osa;
    int len, retval;
    char value[40];
    char addr[80];
    char* msg;

    if (!name  ||  !*name  ||  !*type_name
        ||  !SOCK_HostPortToString(host, port, addr, sizeof(addr))) {
        errno = EINVAL;
        return 0;
    }
    if (!path  ||  !*path)
        path = LBSM_DEFAULT_FEEDFILE;

    if (!(msg = (char*) malloc(20
                               + strlen(name)
                               + strlen(type_name)
                               + strlen(addr)
                               + sizeof(value)))) {
        return 0/*failure*/;
    }
    if (fine)
        NCBI_simple_ftoa(value, min(max(0.0, rate), 100.0), 0);
    else if (LBSM_RERATE_RESERVE <= rate)
        strcpy(value, "RESERVE");
    else if  (rate <= LBSM_RERATE_DEFAULT)
        strcpy(value, "DEFAULT");
    else {
        if (fabs(rate) < SERV_MINIMAL_RATE / 2.0)
            rate = 0.0;
        else if (rate > 0.0)
            rate = max(rate,  SERV_MINIMAL_RATE);
        else if (rate < 0.0)
            rate = min(rate, -SERV_MINIMAL_RATE);
        NCBI_simple_ftoa(value, !rate ? 0.0/*NB: avoid -0.0*/ :
                         min(max(-SERV_MAXIMAL_RATE, rate), SERV_MAXIMAL_RATE),
                         3/*decimal places*/);
    }
    retval = 0;  /*assume worst, a failure*/
    len = sprintf(msg, "%u %s %s%s %s %s\n", (unsigned int) geteuid(),
                  name, fine ? "" : "RERATE ", type_name, addr, value);
    if (len > 0) {
        memset(&sa, 0, sizeof(sa));
        sa.sa_handler = SIG_IGN;
        if (sigaction(SIGPIPE, &sa, &osa) == 0) {
            SOCK cmd;
            SOCK_CreateUNIX(path, 0, &cmd, msg, (size_t)len, fSOCK_LogDefault);
            if (cmd  &&  SOCK_Close(cmd) == eIO_Success)
                retval = 1/*success*/;
            sigaction(SIGPIPE, &osa, 0);
        }
    }
    free(msg);
    return retval;
}
