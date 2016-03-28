#define _GNU_SOURCE
#include "bbapi.h"


pthread_t *bb_proxy_thread;
pthread_cond_t *cond;
pthread_mutex_t *mutex;
pthread_mutex_t *wqtex;
_internal_work_queue_t *wqhead;
_internal_work_queue_t *wqtail;
bool sleeping = false;

static void _mkdir(const char *dir) {
    char *tmp;
    char *p = NULL;
    size_t len;

    tmp = (char*)malloc(strlen(dir)+1);
    strcpy(tmp,dir);

    len = strlen(tmp);
    if(tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            mkdir(tmp, S_IRWXU);
            *p = '/';
        }
}


void spawn_bb_proxy(){
       bb_proxy_thread = (pthread_t *)malloc(sizeof(pthread_t)); 
       cond = (pthread_cond_t *)malloc(sizeof(pthread_cond_t));
       mutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
       wqtex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
       wqhead = NULL;
       wqtail = NULL;

       if (pthread_cond_init(cond,NULL) != 0){
           perror(strerror(errno));
           exit(-1);
       }

       if (pthread_mutex_init(mutex,NULL) != 0){
           perror(strerror(errno));
           exit(-1);
       }
       
       if (pthread_mutex_init(wqtex,NULL) != 0){
           perror(strerror(errno));
           exit(-1);
       }
       pthread_create(bb_proxy_thread,NULL, bb_proxy_func, NULL); 
}


void term_bb_proxy(){
    bool waiting_to_signal = true;

    //No other work items can come after this.
    enqueue_work(1,1,NULL,NULL,NULL,true);

    //We have to make sure this signal is caught by the other thread
    //To avoid an unflushed memfence we check the sleeping variable within
    //the confines of the mutex since lock and unlock forces a memfence.
    while (waiting_to_signal)
    {
        pthread_mutex_lock(mutex);
        if (sleeping){
            waiting_to_signal = false;
        }
        pthread_mutex_unlock(mutex);
    }

    pthread_cond_signal(cond);
    pthread_join(*bb_proxy_thread, NULL);
    pthread_cond_destroy(cond);
    pthread_mutex_destroy(mutex);
    pthread_mutex_destroy(wqtex);
}

void *bb_proxy_func(void *arg){
    bb_bind_cpu();
    bool term_signaled = false;
    struct timespec qtime = {1,0};

    while (!term_signaled){

        _internal_work_queue_t *wi = dequeue_work();

        if (!wi){
            pthread_mutex_lock(mutex);
            sleeping = true;
            pthread_cond_wait(cond,mutex);
            sleeping = false;
            pthread_mutex_unlock(mutex);
        }else{
            if (wi->term_flag){
                term_signaled = true;
            }else{
                char *src = wi->work_item->trans_list->src;
                char *dest = wi->work_item->trans_list->dest;
                move_file(src,dest);
            }

            free(wi);
        }
    }


    return NULL;
}


void move_file(char *src, char *dest){
    sync();
    if (src == NULL || dest == NULL){
        return;
    }
    _mkdir(dest);
    //Theoretically we have a created directory now.
    //Next step -- copy the file.
    int ifd = 0;
    int ofd = 0;
    char buf[4096];
    ssize_t bytes_read = 0;

    if ((ifd=open(src, O_RDONLY) ) == -1){
        fprintf(stderr,"infile %s \n",(strerror(errno)));
        exit(-1);
    }

    if ((ofd=open(dest, O_WRONLY | O_CREAT , 0666)) == -1){
        fprintf(stderr,"outfile %s\n", (strerror(errno)));
        exit(-1);
    }

    while ((bytes_read = read(ifd, buf, sizeof(buf))) > 0){
        char *output = buf;
        ssize_t bytes_written = write(ofd, output, bytes_read);
        if (bytes_written != bytes_read){
            fprintf(stderr,"File move failed\n");
            exit(-1);
        }
    }

    //Successful Termination
    if (bytes_read == 0){
        /* Don't close infile
           It was opened for reading.
           If we close it.. the close call will be intercepted and generate a new
           transfer out of this directory
           
           This is just testing code anyway.
           */
        close(ofd);
    }else{
        perror("Failed to copy file\n");
    }
}

void bb_bind_cpu(){
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(0,&mask);
    sched_setaffinity(0, sizeof(cpu_set_t), &mask);
}


_internal_work_queue_t* dequeue_work(){
    _internal_work_queue_t *ptr = NULL;
    pthread_mutex_lock(wqtex);
    if (wqhead){
        ptr = wqhead;
        wqhead = (_internal_work_queue_t *)wqhead->next;
    }
    pthread_mutex_unlock(wqtex);
    return ptr;
}

int enqueue_work(uint64_t tag, uint32_t num_contrib, uint32_t *contrib, BBTransferDef_t *xfer, BBTransferHandle_t *handle, int term_flag){
    //Slow and ugly... this is in the emulated BBAPI so.. who cares this wont' be in the production code
    pthread_mutex_lock(wqtex);
    if (wqhead == NULL){
        wqhead = (_internal_work_queue_t *)malloc(sizeof(_internal_work_queue_t));
        wqhead->work_item = xfer;
        wqhead->work_handle = handle;
        wqhead->term_flag = term_flag;
        wqhead->next = NULL;        
        wqtail = wqhead;
    }else{
        /* We always grow at the head and pull from the tail
        */
        _internal_work_queue_t *temp = malloc(sizeof(_internal_work_queue_t));
        temp->work_item = xfer;
        temp->work_handle = handle;
        temp->term_flag = term_flag;
        temp->next = NULL;
        wqtail->next = (struct _internal_work_queue_t *)temp;
        wqtail = temp;
    } 
    pthread_mutex_unlock(wqtex);

    pthread_cond_signal(cond);
}

