#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

typedef unsigned long bigInt;
#define LEN 64 // number of bits
#define NMAX 64 // number of elements in array
static const int MAXSIZ=LEN*NMAX;
static const bigInt MSB=0x8000000000000000; //2^64 in hex
const bigInt one=1;
struct LargeNumber{
   bigInt num[NMAX];  //Array of NMAX bigInts
   unsigned int size; //Stores size of number
};
typedef struct LargeNumber LargeNumber;

LargeNumber Zero; // Create a global LN=0
LargeNumber One;  // Create a global LN=1

unsigned int bitCnt(bigInt m);
void zero(LargeNumber *a);
void dec2bin(bigInt n, LargeNumber *a);
unsigned long bin2dec(LargeNumber *a);
void print_bits(LargeNumber* a);
void setbit(LargeNumber *a, unsigned int idx);

LargeNumber copy_number(LargeNumber *a);
LargeNumber add(LargeNumber *a, LargeNumber *b);
LargeNumber subtract(LargeNumber *a, LargeNumber *b);
LargeNumber mult(LargeNumber * a, LargeNumber *b);
LargeNumber square(LargeNumber *a);
LargeNumber divide_2n(LargeNumber *m);
LargeNumber xor(LargeNumber *a, LargeNumber *b);
LargeNumber and(LargeNumber *a, LargeNumber *b);
LargeNumber shift_l(LargeNumber *a, int N, int iFl);
LargeNumber shift_r(LargeNumber *a, int N);
LargeNumber complement(LargeNumber *a);
LargeNumber n_msb_into_new(LargeNumber * r, int n);
LargeNumber n_lsb_into_new(LargeNumber * r, int n);
LargeNumber barrett_reduction(LargeNumber *w, LargeNumber *m);
int test_prime(unsigned long k, unsigned long n, unsigned long aa);
int greater_than(LargeNumber * x, LargeNumber * m);
int equal_to(LargeNumber *a, LargeNumber *b);
int msb(LargeNumber *a);
_Bool isbitset(LargeNumber *a, unsigned long idx);
int msb2(LargeNumber *a, int start);
void prn(LargeNumber *a);




// ======================== Start Main Program ================================

int main() {

// Initialise global LNs
    zero(&Zero);
    dec2bin(1,&One);

    // create number a struct
    LargeNumber a=Zero; // create a LN and initialise to zero
    LargeNumber b=Zero; // create and initialise LN
    LargeNumber r=Zero;
    LargeNumber p=Zero;


    // test divide
    bigInt i;
    bigInt n;
    bigInt f;
    bigInt j;

    // bigInt k = 7;
    // bigInt nn = 4; done

    // bigInt k = 59;
    // bigInt nn = 1085; done

    // bigInt k = 3;
    // bigInt nn = 534; done

    bigInt k = 3;
    bigInt nn = 1103;


    bigInt b_val = 16;
    bigInt a_val = b_val*2;
    // bigInt a_val=3;
    // bigInt y=3;

    a = shift_l(&One, a_val+1, 0);
    a = subtract(&a, &One);
    b = shift_l(&One, b_val, 0);
    b = add(&b, &One);
    printf("%lu, %lu", a_val, b_val);

    // dec2bin(113, &a);
    // dec2bin(114, &b);

    // print_bits(&a);
    // return 0;
    // print_bits(&b);

    clock_t begin, end;
    double time_spent;

    begin = clock();
    barrett_reduction(&a, &b);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f\n", time_spent);

    // for(i=3; i<8; i+=2){
    // printf("a value %lu\n", i);
    //    if(test_prime(k, nn, i)){
    //         printf("yes\n" );
    //         exit(0);
    //     }
    // }


    // print_bits(&r);
    // printf("r is %d\n", msb(&r));


return 0;
}

// ======================== End Main Program ================================

int test_prime(unsigned long k, unsigned long n, unsigned long aa){
    LargeNumber r=Zero;
    LargeNumber a=Zero;
    LargeNumber p = Zero;
    LargeNumber k_ = Zero;
    LargeNumber twon = One;
    double time_spent;

    dec2bin(k, &k_);
    dec2bin(aa, &a);
    twon = shift_l(&One, n, 1);

    clock_t begin = clock();

    p = mult(&k_, &twon);
    p = add(&p, &One);

    clock_t end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("make p %f\n", time_spent);

    begin = clock();
    // find a^k
    r = square(&a);
    bigInt i;
    for(i=0; i<k-2; i++){
        r = mult(&r, &a);
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("squaring %f\n", time_spent);


    // frist barrett_reduction
    begin = clock();
    r = barrett_reduction(&r, &p);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("first br %f\n", time_spent);

    // square and barrett_reduction n-1 times
    for(i=0; i<n-1; i++){
        begin = clock();
        r = square(&r);
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("loop squaring %f\n", time_spent);
        begin = clock();
        r = barrett_reduction(&r, &p);
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("this is the %lu step %f\n", i, time_spent);


    }

    p = subtract(&p, &One);

    // print_bits(&p);
    // printf("%lu\n", bin2dec(&p));
    // print_bits(&r);
    // printf("result %lu\n", bin2dec(&r));

    if(equal_to(&r, &p)){
        return 1;
    }
    else{
        return 0;
    }
}

void zero(LargeNumber *a){
  for(int i=0;i<NMAX;i++){
    a->num[i]=0;
  }
  a->size=0;
  return;
}

void dec2bin(bigInt n, LargeNumber *a){
  a->num[0]=n;
  a->size=bitCnt(n);
  return;
}

unsigned int bitCnt(bigInt m){
  unsigned int n=LEN;
  while((m & MSB)==0){
    n--;
    m=m<<1;
  }
  return n;
}

void print_bits(LargeNumber* a){
  int i;
  unsigned int j,jMax;
  if (a->size == 0){
    jMax=0;
  }else{
    jMax=(NMAX-1)/LEN;
  }
  for (j=0; j < NMAX; j++){
    i = LEN;
    while (i != 0){
      i--;
      printf("%lu", (a->num[j] >> i) & 0x01);
    }
    printf("\n");
  }
}

unsigned long bin2dec(LargeNumber *a){
  bigInt m=0;
  if(a->size<=LEN){
    m=a->num[0];
  }else{
    printf("bin2dec: number too big to print\n");
  }
  return m;
}

LargeNumber copy_number(LargeNumber *a){
  LargeNumber r=*a;
  return r;
}

LargeNumber shift_l(LargeNumber *a, int N, int iFl){
    long i, s, j, carry;
    // check that bitshift is possible given fixed max size of numbers
    if (iFl){
      int ok=(a->size+N)<=MAXSIZ;
      if (!ok){
        printf("aborting shift_l output too long %d\n", N);
        exit(0);
      }
    }
    // copy a to new LN r
    LargeNumber r=*a;


    for(s=0; s<N; s++){ // do it N times
      // make shift by 1, storing carry bits temporarily as necessary
      carry=0;
      for(j=0; j<NMAX; j++){
          if (r.num[j]&MSB) carry+=(one<<j); // each bit of carry holds the carry bit for one member of num array.
          r.num[j] = r.num[j] << 1;
      }

      // if there was a one at end put it in next element
      for(i=1; i<NMAX; i++){ // this is only for NMAX-1 members of array num, since we can't extend the array and we start on number 1 (not 0). The MS carry bit is lost
          if(carry&one){
            r.num[i]=r.num[i]|1;//sets LSB of ith element of num
          }
          carry=carry>>1;// move along one in carry for next iteration
      }
    }
    // don't forget to increment size of r
    r.size=r.size+N;
    return r;
}

LargeNumber shift_r(LargeNumber *a, int N){
    long i, s, j, carry;
    // copy a to new LN r
    LargeNumber r=*a;

    for(s=0; s<N; s++){ // do it N times
      // make shift by 1, storing carry bits temporarily as necessary
      carry=0;
      for(j=0; j<NMAX; j++){
          if (r.num[j]&1) carry+=(one<<j); // each bit of carry holds the carry bit for one member of num array.
          r.num[j] = r.num[j] >> 1;
      }
      // if there was a one at end put it in next element
      for(i=0; i<NMAX-1; i++){ // this is only for NMAX-1 members of array num, since we can't extend the array and we start on number 1 (not 0). The LS carry bit is lost
          carry=carry>>1;// move along one in carry for next iteration (lose the LS one)
          if(carry&1){
            // printf("%d\n",one_storage[i]);
            //setbit(r->num, one_storage[i]+1);
            r.num[i]=r.num[i]|MSB;//sets MSB of ith element of num
          }
      }
    }
    // don't forget to increment size of r
    r.size=r.size-N;
    return r;
}

// sets individual bit to 1
void setbit(LargeNumber *a, unsigned int idx){
  unsigned int i, iBit, idy;
  //field[idx / 8] |= 1u << (idx % 8);
  i=idx/LEN; // which num is it in?
  iBit=(idx % LEN); // assumes idx is numbered from 0
  a->num[i]=a->num[i] |= 1 << iBit;
  idy=idx+1;
  if(idy > a->size) a->size=idy;
}
// void setbit(LargeNumber *a, unsigned int idx){
//   a->num[idx / LEN] |= 1 << (idx % LEN);
// }
/*
// sets individual bit to 0
void unsetbit(unsigned long *field, unsigned int idx) {
  field[idx / 8] &= ~(1u << (idx % 8));
}

// checks if individual bit is 1
_Bool isbitset(unsigned long *field, unsigned long idx) {
  return field[idx / 8] & (1u << (idx % 8));
}
*/

// checks if individual bit is 1
_Bool isbitset(LargeNumber *a, unsigned long idx) {
  unsigned int i, iBit;
  i=idx/LEN; // which num is it in?
  iBit=(idx % LEN); // assumes idx is numbered from 0
  return (a->num[i] & (one << iBit));
}
/*

*/
int is_number_zero(LargeNumber *a){
// Like msb, this should only be used if you don't have a length for the number.
// If the number is zero, it's length should be zero, and vice versa.
  unsigned int i;
  for(i=0; i<=a->size/LEN; i++){
    if(a->num[i] != 0){
      return 0;
    }
  }
  return 1; //meaning number is zero
}

int msb(LargeNumber *a){
// A function to find the msb when we don't know the size of the LargeNumber
// If we do know the size, then msb is just size-1 which is faster than this function
  int i, j, max;
  max = MAXSIZ; // starting value
  for (j=NMAX-1; j>=0; j--){ // looping through elements of array
    for(i=0;i<LEN;i++){ // looping through bits in element starting at MSB
      if((a->num[j] << i) & MSB){
        return max-i-1;
      }
    }
    max=max-LEN;
  }
  return -1;
}

LargeNumber add(LargeNumber *a, LargeNumber *b){
    LargeNumber carry  = and(a, b);
    LargeNumber result = xor(a, b);
    LargeNumber shifted_carry=Zero;

    while(carry.size != 0) // while number is not zero
    {
        shifted_carry = shift_l(&carry, 1, 0);
        carry = and(&result, &shifted_carry);
        result = xor(&result, &shifted_carry);
    }
    return result;
}

LargeNumber divide_2n(LargeNumber *m){
    bigInt i,n;
    n = m->size;
    LargeNumber mp = shift_l(m, n, 1);
    LargeNumber R = One;
    LargeNumber r = shift_l(&One, 2*n, 1);
    r = subtract(&r,&mp);


    for(i=0; i<n; i++){

        r = shift_l(&r, 1, 1);

        R = shift_l(&R, 1, 1);

        if(greater_than(&r, &mp) || equal_to(&r, &mp)){ //shouldn't this be >=?
            // print_bits(r);


            r = subtract(&r, &mp);
            R = add(&R, &One);
        }

    }
    return R;
}

LargeNumber subtract(LargeNumber *a, LargeNumber *b){
    // We are going to do this a different way, namely using
    // bitwise complement and add
    // Be careful - works only for b<a, but must be true for subtract in unsigned arithmetic
    int ok=greater_than(a,b);
    if (!ok){
      printf("aborting subtract since difference negative\n");
      exit(0);
    }
//    int larger = ((a->size) >= (b->size)) ? (a->size) : (b->size); //larger one
    LargeNumber r=complement(b);
    // printf("ayyyy %lu \n \n", bin2dec(b));
    LargeNumber bp=add(&r,&One); // this should be "-b" according to rules of twos complement
    r=add(a,&bp);

    r.size=msb(&r)+1;

    return r;
}

LargeNumber mult(LargeNumber * a, LargeNumber *b){
    // if a or b equal to zero return Zero
    if(equal_to(a, &Zero) || equal_to(a, &Zero)){
        return Zero;
    }
    unsigned int i, j;
    LargeNumber aShft = *a;

    LargeNumber bShft = *b;
    LargeNumber result = Zero;
    // printf("size is %d\n", a->size+b->size);
    int ok=(a->size+b->size <= MAXSIZ);
    if (!ok){
        // prn(b);
      printf("aborting mult since output too long\n");
      exit(0);
    }
    j=0;
    for(i=0; i < a->size ; i++){ // loop over bits in a
        // printf("a size %lu\n", a->size);
        if(aShft.num[j]&1){
            result=add(&result,&bShft);
        };
        aShft.num[j]=aShft.num[j]>>1; // just shift the current element of a
        bShft=shift_l(&bShft,1,1); // shift the whole of b
        if( (i+1) % LEN == 0 ) j++; // change which element of num array in a we are working with
    }
    return result;
}

LargeNumber square(LargeNumber * a){
    LargeNumber result = mult(a, a);
    return result;
}

LargeNumber and(LargeNumber *a, LargeNumber *b){
    int i;
    int smaller = ((a->size) <= (b->size)) ? (a->size) : (b->size); //smaller one
    //int larger  = ((a->size) >= (b->size)) ? (a->size) : (b->size); //larger one
    LargeNumber result=Zero;
    for(i=0; i<=smaller/LEN; i++){
      result.num[i] = a->num[i] & b->num[i];
    }
    result.size=msb(&result)+1;
    return result;
}

LargeNumber xor(LargeNumber *a, LargeNumber *b){
    int i;
    int larger = ((a->size) >= (b->size)) ? (a->size) : (b->size); //larger one
    LargeNumber result=Zero;
    result.size = larger; //size of bigger one
    for(i=0; i<=(larger-1)/LEN; i++){
      result.num[i] = a->num[i] ^ b->num[i];
    }
    // printf("%d\n", result.size );
    // if ((a->size) == (b->size)){
    //     result.size=msb(&result)+1;
    //     printf("%d\n", result.size );
    // }
    return result;
}

LargeNumber complement(LargeNumber *a){
    unsigned int i;
    LargeNumber r = *a;
    for(i=0; i < NMAX; i++){
        r.num[i] = ~(r.num[i]);
    }
    if(a->size < MAXSIZ ){
      r.size=MAXSIZ;
    }else{
      r.size=msb(&r)+1;
    }
    return r;
}

LargeNumber n_msb_into_new(LargeNumber * r, int n){
    int rSiz, diff;
    LargeNumber r_copy = *r;
    rSiz=r_copy.size;
    diff=rSiz-n;
    r_copy = shift_r(&r_copy, diff);
    return r_copy;
}

LargeNumber n_lsb_into_new(LargeNumber * r, int n){
    int i;

    LargeNumber mask = One;
    for (i=1; i<n; i++){
      mask=shift_l(&mask,1,1);
      mask=xor(&mask,&One);
    }
    LargeNumber result = and(r,&mask);

    return result;
}


int equal_to(LargeNumber *a, LargeNumber *b){
    int i;
    for(i=0; i<NMAX; i++){
        if(a->num[i] != b->num[i]){
            return 0;
        }
    }
    return 1;
}


int greater_than(LargeNumber * x, LargeNumber * m){
    LargeNumber xor_xm = xor(x, m);
    unsigned long msb_1 = msb(&xor_xm);
    if(isbitset(x, msb_1)){
        return 1;
    }
    else{
        return 0;
    }
}

LargeNumber barrett_reduction(LargeNumber *w, LargeNumber *m){
    int n;
    n=m->size;
    clock_t begin;
    clock_t end;
    double time_spent;

    begin = clock();
    LargeNumber r = divide_2n(m); // 2^(2n)/m
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f", time_spent);

    begin = clock();
    LargeNumber y = mult(w, &r); // multiply w and r into y
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f", time_spent);

    begin = clock();
    LargeNumber q = shift_r(&y, 2*n);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f", time_spent);

    begin = clock();
    LargeNumber z = mult(&q, m);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f", time_spent);

    begin = clock();
    LargeNumber x = subtract(w,&z); // subtract
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(", %f", time_spent);

    if(greater_than(&x,m) || equal_to(&x,m)){
      begin = clock();
      x=subtract(&x,m);
      end = clock();
      time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      printf("\tbarrett subtract loop %f\n", time_spent);
    }
    return x;
}

void prn(LargeNumber *a){
    print_bits(a);
    printf("Size = %d\n",a->size);
    printf("%lu\n", bin2dec(a));
    printf("\n");
}
