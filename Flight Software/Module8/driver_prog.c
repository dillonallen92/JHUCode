#include <stdio.h>
#include <stdlib.h>

typedef struct mp_t{
    int userData;
};

struct mp_t *MP_Create(int blocksize, int num_blocks);
// int MP_AllocBlock();
// int MP_FreeBlock(int block_idx);
// void *MP_GetBlockData(int block_idx);
// int MP_Destroy(mp_t **mp);

int main()
{
    int mySize = 10;
    struct mp_t *myMemPoolPtr = MP_Create(sizeof(int), mySize);
    for(int i = 0; i < mySize; i++)
    {
        myMemPoolPtr[i].userData = i;
    }
    for(int i = 0; i < mySize; i++)
    {
        printf("%d \n", myMemPoolPtr[i].userData);
    }

}

struct mp_t *MP_Create(int blocksize, int num_blocks)
{
    if(blocksize > 1024)
    {
        printf("Max blocksize is 1024");
        return NULL;
    }

    if(num_blocks > 100)
    {
        printf("Max number of blocks is 100");
        return NULL;
    }

    struct mp_t *myMP = (struct mp_t*) malloc(num_blocks * (size_t) blocksize);

    return myMP;
}