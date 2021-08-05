#include <stdio.h>
#include <stdlib.h>

typedef struct mp_t{
    short int isFree;
    short int userData;
} myMP;

struct mp_t *MP_Create(int blocksize, int num_blocks);
int MP_AllocBlock(struct mp_t *myMP, int num_blocks);
int MP_FreeBlock(struct mp_t *myMP, int block_idx);
void *MP_GetBlockData(struct mp_t *myMp, int block_idx);
int MP_Destroy(struct mp_t **myMP);

int main()
{
    int mySize = 100;
    struct mp_t *myMemPoolPtr = MP_Create(2, mySize);

    printf("Allocating %d blocks with %d bytes: \n", mySize, 2);
    if(myMemPoolPtr == NULL)
    {
        printf("Improper Size");
    }else{
        for(int i = 0; i < mySize; i++)
        {
        printf("%d \n", myMemPoolPtr[i].userData);
        }
    }

    printf("Attempting to take block 101.. Freeing first\n");
    if(MP_FreeBlock(myMemPoolPtr, 101))
    {
        printf("Block 101 freed successfully! \n");
    }else
    {
        printf("Error freeing block. Out of bounds. \n");
    }

    printf("Freeing memory...\n");
    
    if(MP_Destroy(myMemPoolPtr))
        printf("memory cleared! \n");
    else
        printf("Error clearning memory;");
    

    /*

    ** THIS IS CODE TO MAKE SURE EVERYTHING WORKED **

    printf("Index of next allocated block: %d \n", MP_AllocBlock(myMemPoolPtr, mySize));
    printf("Freeing block 3\n");
    if(MP_FreeBlock(myMemPoolPtr, 2))
    {
        printf("Successfully freed block 3\n");
    }else{
        printf("Error in freeing block 3\n");
    }
    printf("Index of next freed block: %d \n", MP_AllocBlock(myMemPoolPtr, mySize));
    printf("User data in block 4: ");
    if(MP_GetBlockData(myMemPoolPtr, 3) != NULL)
    {
        printf("%d", (int*)MP_GetBlockData(myMemPoolPtr, 3));
    }else{
        printf("\nFree block, no user data.\n");
    }

    printf("\nFreeing memory...\n");
    
    if(MP_Destroy(myMemPoolPtr))
        printf("memory cleared! \n");
    else
        printf("Error clearning memory;");
    */

    return 0;
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

    struct mp_t *myMP = (struct mp_t*) calloc(num_blocks, (size_t) blocksize);
    
    for(int i = 0; i < num_blocks; i++)
    {
        myMP[i].isFree = 0;
        myMP[i].userData = i+1;
    }

    return myMP;
}

int MP_AllocBlock(struct mp_t *myMP, int num_blocks)
{
    int idx = 0;
    while(myMP[idx].isFree != 1 && idx < num_blocks)
    {
        idx++;
    }
    if(idx >= num_blocks)
    {
        idx = -1;
    }else{
        myMP[idx].isFree = 0;
    }
    return idx;
}

int MP_FreeBlock(struct mp_t *myMP, int block_idx)
{
    int didSucceed = 0;
    if(myMP[block_idx].isFree == 1)
    {
        printf("Block is already free.");
    }else if(myMP[block_idx].isFree == NULL)
    {
        didSucceed = 0;
    }else{
        myMP[block_idx].userData = 0;
        myMP[block_idx].isFree = 1;
        didSucceed = 1;
    }

    return didSucceed;
}

void *MP_GetBlockData(struct mp_t *myMP, int block_idx)
{
    if(myMP[block_idx].isFree == 0)
        return (void*) (myMP + block_idx)->userData;
    else
        return NULL;
}

int MP_Destroy(struct mp_t **myMP)
{
    free(myMP);
    myMP = NULL;

    // check to see if it actually was set to null or if error
    if(myMP == NULL)
    {
        return 1;
    }else
    {
        return 0;
    }
}