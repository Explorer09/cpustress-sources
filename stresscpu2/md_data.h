#ifndef _md_data_h_
#define _md_data_h_

typedef struct
{
    int       nri;
    int *     iinr;
    int *     gid;
    int *     shift;
    int *     jindex;
    int *     jjnr;
    int       count;
} 
nlist_t;


typedef struct
{
    int            natoms;
    float *        x;
    float *        shiftvec;
    float *        charge;
    int *          type;
    float *        nbparam;
    int            ntype;
    nlist_t        nlist;
} 
md_data_t;


int
setup_md_data(md_data_t *     data);



void
release_md_data(md_data_t *   data);



#endif

