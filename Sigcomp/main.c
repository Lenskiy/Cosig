//
//  main.c
//  Sigcomp
//
//  Created by Artem Lenskiy on 15/03/2016.
//  Copyright © 2016 Artem Lenskiy. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846264338327950288
#endif

#define BLOCK_SIZE 16

static long current_sample_in = 0;
static long current_sample_out = 0;
static long total_n_samples = 0;

//===============================================================================================================
long readTextFile(char path[], short *samples){
    
    FILE *fd = NULL;
    fd = fopen(path, "r");
    
    if(fd == NULL){         //the file could be opened.
        return -1;
    }
    
    while(fscanf(fd, "%hd", samples + total_n_samples++) != EOF);
    total_n_samples--;
    
    fclose(fd);
    
    return total_n_samples;}

//---------------------------------------------------------------------------------------------------------------
long writeTextFile(char path[], const short *samples, unsigned long length){
    
    FILE *fd = NULL;
    fd = fopen(path, "w+");
    
    if(fd == NULL){         //the file could be opened.
        return -1;
    }
    
    unsigned long i = 0;
    while(i < length);
    fprintf(fd, "%hd", *(samples + i++));
    
    fclose(fd);
    
    return total_n_samples;}

//---------------------------------------------------------------------------------------------------------------
unsigned getBlock(short *samples, short *block){
    
    unsigned n_samples_read = 0;
    
    while (n_samples_read < BLOCK_SIZE && current_sample_in < total_n_samples) // replace with memcpy
        block[n_samples_read++] = samples[current_sample_in++];

    return n_samples_read;}

//---------------------------------------------------------------------------------------------------------------
unsigned putBlock(short *samples, const short *block){
    
    unsigned n_samples_written = 0;
    
    while (n_samples_written < BLOCK_SIZE && current_sample_out < total_n_samples) // replace with memcpy
        samples[current_sample_out++] = block[n_samples_written++];
    
    return n_samples_written;}

//===============================================================================================================
void printBlock(short *block){
    unsigned int i = 0;
    for(;i < BLOCK_SIZE; i++)
        printf("%hd, ", block[i]);
    printf("\n");
    
    return;}

//===============================================================================================================
void quantizeBlock(short *block){
    unsigned int i = 0;
    for(;i < BLOCK_SIZE; i++)
        block[i] = block[i] >> 1;
    
    return;}
//---------------------------------------------------------------------------------------------------------------
void invQuantizeBlock(short *block){
    unsigned int i = 0;
    for(;i < BLOCK_SIZE; i++)
        block[i] = block[i] << 1;
    
    return;}
//===============================================================================================================
//input: f, N; output: F
short dct_direct(double *f, double *F ){
    double a[BLOCK_SIZE], sum;
    short i, v;

    a[0] =  sqrt(1.0 / BLOCK_SIZE); // sqrt ( 1.0 / N );
    for ( i = 1; i < BLOCK_SIZE; ++i )
        a[i] =   sqrt(2.0 / BLOCK_SIZE) ; // sqrt ( 2.0 / N );
    
    for ( v = 0; v < BLOCK_SIZE; ++v ) {
        sum = 0.0;
        for ( i = 0; i < BLOCK_SIZE; ++i ) {
            sum += f[i] * cos( (2*i+1) * v * M_PI/(2*BLOCK_SIZE) );
        } //for j
        F[v] = a[v] * sum;
    } //for v
    return 1;}
//---------------------------------------------------------------------------------------------------------------
//f[i][j] * coef
//input: N, F; output f
short idct_direct(double *F, double *f ){
    double a[BLOCK_SIZE], sum;
    short i, v;

    a[0] =  sqrt(1.0 / BLOCK_SIZE);
    for ( i = 1; i < BLOCK_SIZE; ++i )
        a[i] =  sqrt(2.0 / BLOCK_SIZE);
    
    for ( i = 0; i < BLOCK_SIZE; ++i ) {
        sum = 0.0;
        for ( v = 0; v < BLOCK_SIZE; ++v ) {
            sum += a[v]*F[v]*cos( (2*i+1) * v * M_PI/(2*BLOCK_SIZE) );
        } //for v
        f[i] =  sum;
    } //for i
    return 1; }

//---------------------------------------------------------------------------------------------------------------
// change values from short to double and vice versa.
short dct (const short *f, short *F ){
    double  temp_f[BLOCK_SIZE] = {0}, temp_F[BLOCK_SIZE] = {0};
    int  i;

    for ( i = 0; i < BLOCK_SIZE; ++i )
        temp_f[i] = (double) f[i];
    
    dct_direct (temp_f, temp_F); //DCT operation
    
    for ( i = 0; i < BLOCK_SIZE; ++i )
        F[i] = (short ) ( floor (temp_F[i] + 0.5) ); //rounding
    
    return 1;}
//---------------------------------------------------------------------------------------------------------------
// change values from short to doulbe, and vice versa.
short idct (const short *F, short *f ){
    double  temp_f[BLOCK_SIZE] = {0}, temp_F[BLOCK_SIZE] = {0};
    int  i;

    for ( i = 0; i < BLOCK_SIZE; ++i )
        temp_F[i] = (double) F[i];
    
    idct_direct (temp_F, temp_f );        //IDCT operation
    
    for ( i = 0; i < BLOCK_SIZE; ++i )
        f[i] = (short ) ( floor (temp_f[i] + 0.5) ); //rounding
    
    return 1; }
//===============================================================================================================
short sub (const short a, const short b){
    return (a - b);
}
//---------------------------------------------------------------------------------------------------------------
void applyOpBlock(short *block1, const short *block2, short (*op) (short, short)){
    unsigned i = 0;
    for(;i < BLOCK_SIZE; ++i)
        block1[i] = op(block1[i], block2[i]);
}
//===============================================================================================================
int main(int argc, const char * argv[]) {
    
    char path[] = "/Users/artemlenskiy/Documents/Research/Matlab/HRVformHWsensor/finger pulse.txt";
    
    short *samples = (short *)malloc(1e6 * sizeof(short)); //Allocted memory to store 1 million of samples of size short
    
    if(readTextFile(path, samples) < 0){
        puts("Couldn't read file");
        return -1;
    }
    
    long n_samplesInBlock = 0;
    short *block =      (short *) malloc(BLOCK_SIZE * sizeof(short));
    short *tempBlock =  (short *) malloc(BLOCK_SIZE * sizeof(short));
    short *dctBlock =   (short *) malloc(BLOCK_SIZE * sizeof(short));
    short *reconBlock = (short *) malloc(BLOCK_SIZE * sizeof(short));
    
    //test coding and decoding
    do{
        n_samplesInBlock = getBlock(samples, block);
        memcpy(tempBlock, block, BLOCK_SIZE * sizeof(short));
        //printBlock(tempBlock);
        quantizeBlock(tempBlock);
        dct(tempBlock, dctBlock);
        //printBlock(dctBlock);
        idct(dctBlock, reconBlock);
        invQuantizeBlock(reconBlock);
        applyOpBlock(block, reconBlock, sub);
        printBlock(block);
    }while(n_samplesInBlock != 0);
    
    free(reconBlock);
    free(dctBlock);
    free(tempBlock);
    free(block);
    free (samples);
    
    return 0;
}
