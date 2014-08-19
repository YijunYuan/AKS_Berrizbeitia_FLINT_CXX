/** Implementation of AKS-Berrizbeitia Algorithm
 *  Coder: Yijun Yuan
 *  Email: 941201yuan@gmail.com
 *  This file depends on the William Hart's FLINT
 *  The implementation of original AKS Algorithm made by me can be found at here:
 *      https://github.com/YijunYuan/AKS_FLINT
 *
*/
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fmpz_mod_poly.h>
#include<flint/arith.h>
#include"IQF_CALC.h"
#include<math.h>

const bool PRIME_TABLE[] ={ 0,
0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
bool CASE_I(fmpz_t n){
    /***********************     Generate a   **************************/
    mpz_t e;    mpz_init(e);fmpz_get_mpz(e,n);//e=n
    mpz_t a_mpz;mpz_init(a_mpz);
    for(mpz_set_ui(a_mpz,2);;mpz_add_ui(a_mpz,a_mpz,1)){
        if(mpz_jacobi(a_mpz,e)==-1)//Generate a
            break;
    }
    fmpz_t a;    fmpz_init(a);fmpz_set_mpz(a,a_mpz);
    mpz_clear(a_mpz);
    /**************************Generate s, k****************************/
    ulong temp_ui;
    fmpz_t A;    fmpz_init_set(A,n);
    fmpz_sub_ui(A,A,1);//A=n-1
    const ulong k=fmpz_val2(A);//k=v2(n-1)
    mpfr_t n_mpfr;mpfr_init2(n_mpfr,200);
        mpfr_t t_mpfr;mpfr_init2(t_mpfr,200);
    mpfr_set_z(n_mpfr,e,MPFR_RNDD);//n=e
    mpz_clear(e);
    mpfr_log2(n_mpfr,n_mpfr,MPFR_RNDD);
        mpfr_div_ui(t_mpfr,n_mpfr,2,MPFR_RNDD);
        temp_ui=mpfr_get_ui(t_mpfr,MPFR_RNDU);//temp_ui=ceil((1/2)log(n))
        mpfr_clear(t_mpfr);
    mpfr_log2(n_mpfr,n_mpfr,MPFR_RNDD);
    mpfr_mul_ui(n_mpfr,n_mpfr,2,MPFR_RNDD);
    const ulong s=mpfr_get_ui(n_mpfr,MPFR_RNDU);//s=ceil(2log(log(n)))
    mpfr_clear(n_mpfr);
    /******************************Step 1*******************************/
    fmpz_fdiv_q_ui(A,A,2);//A=(n-1)/2
    fmpz_powm(A,a,A,n);//A=a^((n-1)/2) mod n
    fmpz_add_ui(A,A,1);
    if(fmpz_equal(A,n)==0){
        fmpz_clear(a);fmpz_clear(A);
        return false;
    }
    if(k>=temp_ui){
        fmpz_clear(a);fmpz_clear(A);
        return true;
    }
    /******************************Step 3*******************************/
    fmpz_t temp;    fmpz_init_set_ui(temp,1);
    fmpz_t m;       fmpz_init_set_ui(m,1);
    ulong card_max,card1=0,card2=0,i;
    if(s<k)
        card_max=0;
    else
        card_max=s-k;
    card_max=exp2l(card_max);
    fmpz_t* Set1=NULL;fmpz_t* Set2=NULL;//Generate set1 and set2
    Set1=(fmpz_t*)flint_calloc(card_max,sizeof(fmpz_t));//|S |=2^max(s-k,0)
    Set2=(fmpz_t*)flint_calloc(card_max,sizeof(fmpz_t));//|S'|=2^max(s-k,0)
    fmpz_set(Set1[0],m);fmpz_set(Set2[0],m);
    card1++;            card2++;
    bool finded;
    fmpz_t exp2_k;fmpz_init_set_ui(exp2_k,2);
    fmpz_pow_ui(exp2_k,exp2_k,k);//exp2_k=2^k
    while(card1<card_max){
         while(1){
            finded=false;
            fmpz_powm(A,m,exp2_k,n);//A=m^exp2_k mod n
            for(i=0;i<card2;i++){
                if(fmpz_equal(A,Set2[i])==1){
                    fmpz_add_ui(m,m,1);finded=true;break;
                }
            }
            if(finded==false)break;
        }
        fmpz_mul_ui(temp,exp2_k,card1);
        fmpz_add_ui(temp,temp,1);
        if(fmpz_cmp(m,temp)>0){
            fmpz_clear(a);
            fmpz_clear(A);
            fmpz_clear(temp);
            fmpz_clear(m);
            fmpz_clear(exp2_k);
            flint_free(Set1);flint_free(Set2);
            return false;
        }
        fmpz_gcd(temp,m,n);
        if(fmpz_cmp_ui(temp,1)>0){
            fmpz_clear(a);
            fmpz_clear(A);
            fmpz_clear(temp);
            fmpz_clear(m);
            fmpz_clear(exp2_k);
            flint_free(Set1);flint_free(Set2);
            return false;
        }
        for(i=0;i<card2;i++){
            fmpz_sub(temp,A,Set2[i]);
            fmpz_gcd(temp,temp,n);
            if(fmpz_cmp_ui(temp,1)>0){
                fmpz_clear(a);
                fmpz_clear(A);
                fmpz_clear(temp);
                fmpz_clear(m);
                fmpz_clear(exp2_k);
                flint_free(Set1);flint_free(Set2);
                return false;
            }
        }
        for(i=0;i<card1;i++){
            if(fmpz_equal(Set1[i],m)==1){
                goto ext1;
            }
        }
        fmpz_set(Set1[card1],m);card1++;
        ext1:;
        for(i=0;i<card2;i++){
            if(fmpz_equal(Set2[i],A)==1){
                goto ext2;
            }
        }
        fmpz_set(Set2[card2],A);card2++;
        ext2:;
    }
    flint_free(Set2);
    fmpz_clear(exp2_k);
    /*******************************Step 4*******************************/
    fmpz_fdiv_q_2exp(temp,n,s);//temp=[n/2^s]
    fmpz_powm(A,a,temp,n);//A=a^[n/2^s] mod n
    fmpz_fdiv_r_2exp(temp,n,s);//temp=n%2^s
    temp_ui=fmpz_get_ui(temp);//temp_ui=n%2^s
    ulong exp2_s=exp2l(s);//s=2^s
    fmpz_mul_si(a,a,-1);
    fmpz_mod_poly_t modulo;fmpz_mod_poly_init(modulo,n);
    fmpz_mod_poly_set_coeff_fmpz(modulo,0,a);
    fmpz_mod_poly_set_coeff_ui  (modulo,exp2_s,1);//modulo=x^2^s-a
    fmpz_mod_poly_t poly1;fmpz_mod_poly_init(poly1,n);
    fmpz_mod_poly_set_coeff_ui(poly1,0,1);
    fmpz_mod_poly_t poly2;fmpz_mod_poly_init(poly2,n);
    fmpz_mod_poly_set_coeff_ui(poly2,0,1);
    fmpz_mod_poly_t polytmp;fmpz_mod_poly_init(polytmp,n);
    for(i=0;i<card1;i++){
        fmpz_mod_poly_set_coeff_fmpz(poly1,1,Set1[i]);//poly1=1+m*x
        fmpz_mod_poly_powmod_fmpz_binexp(polytmp,poly1,n,modulo);//e=poly1^n MOD modulo
        fmpz_mul(temp,Set1[i],A);//temp=m*a^[n/2^s] mod n
        fmpz_mod_poly_set_coeff_fmpz(poly2,temp_ui,temp);//poly2=1+m*a^[n/2^s]*x^(n%2^s)
        if(fmpz_mod_poly_equal(polytmp,poly2)==0){
            fmpz_clear(a);
            fmpz_clear(A);
            fmpz_clear(temp);
            fmpz_clear(m);
            fmpz_mod_poly_clear(poly1);
            fmpz_mod_poly_clear(poly2);
            fmpz_mod_poly_clear(polytmp);
            flint_free(Set1);
            return false;
        }
    }
    fmpz_clear(a);
    fmpz_clear(A);
    fmpz_clear(temp);
    fmpz_clear(m);
    fmpz_mod_poly_clear(poly1);
    fmpz_mod_poly_clear(poly2);
    fmpz_mod_poly_clear(polytmp);
    flint_free(Set1);
    return true;
}
bool CASE_II(fmpz_t n){

    /***********************     Generate a   **************************/
    mpz_t e;    mpz_init(e);fmpz_get_mpz(e,n);//e=n
    mpz_t a_mpz;mpz_init(a_mpz);
    mpz_t temp_mpz;mpz_init(temp_mpz);
    for(mpz_set_ui(a_mpz,2);;mpz_add_ui(a_mpz,a_mpz,1)){
        if(mpz_jacobi(a_mpz,e)==-1){//Generate a
            mpz_set_ui(temp_mpz,1);
            mpz_sub(temp_mpz,temp_mpz,a_mpz);
            if(mpz_jacobi(temp_mpz,e)==-1){
                mpz_clear(temp_mpz);
                break;
            }
        }
    }
    fmpz_t a;    fmpz_init(a);fmpz_set_mpz(a,a_mpz);
    mpz_clear(a_mpz);
    /**************************Generate k,t****************************/
    ulong temp_ui;
    fmpz_t A;    fmpz_init_set(A,n);
    fmpz_add_ui(A,A,1);//A=n+1
    const ulong k=fmpz_val2(A);//k=v2(n+1)
    mpfr_t n_mpfr;mpfr_init2(n_mpfr,256);
        mpfr_t t_mpfr;mpfr_init2(t_mpfr,256);
    mpfr_set_z(n_mpfr,e,MPFR_RNDD);//n=e
    mpz_clear(e);
    mpfr_log2(n_mpfr,n_mpfr,MPFR_RNDD);
        mpfr_div_ui(t_mpfr,n_mpfr,2,MPFR_RNDD);
        temp_ui=mpfr_get_ui(t_mpfr,MPFR_RNDU);//temp_ui=ceil((1/2)log(n))
        mpfr_clear(t_mpfr);
    mpfr_log2(n_mpfr,n_mpfr,MPFR_RNDD);
    mpfr_mul_ui(n_mpfr,n_mpfr,2,MPFR_RNDD);
    const ulong t=mpfr_get_ui(n_mpfr,MPFR_RNDU)+1;//t=ceil(2log(log(n)))+1
    mpfr_clear(n_mpfr);
    /******************************Step 1*******************************/
    /** (a) **/
    fmpz_sub_ui(A,n,1);//A=n-1
    fmpz_fdiv_q_ui(A,A,2);//A=(n-1)/2
    fmpz_powm(A,a,A,n);//A=a^((n-1)/2) mod n
    fmpz_add_ui(A,A,1);
    if(fmpz_equal(A,n)==0){
        fmpz_clear(a);fmpz_clear(A);
        return false;
    }
    /** (b) **/
    fmpz_mul_si(A,a,-1);
    fmpz_add_ui(A,A, 1);//A=1-a
    fmpz_t uin,vin,uout,vout;
    fmpz_init_set_ui(uin,1);
    fmpz_init_set_ui(vin,1);
    fmpz_init(uout);
    fmpz_init(vout);
    IQF_EXP(n,A,n,uin,vin,uout,vout);///ERROR in this function!!!
    fmpz_clear(uin);fmpz_clear(vin);fmpz_sub(vout,vout,n);
    if(fmpz_equal_si(uout,1)==0||fmpz_equal_si(vout,-1)==0){
        fmpz_clear(uout);fmpz_clear(vout);
        fmpz_clear(a);fmpz_clear(A);
        return false;
    }
    fmpz_clear(uout);fmpz_clear(vout);
    /** (c) **/
    if(k>=temp_ui){
        fmpz_clear(a);fmpz_clear(A);
        return true;
    }

    /******************************Step 3*******************************/
    ulong card_max;
    if(t>k)
        card_max=t-k;
    else
        card_max=0;
    card_max=exp2l(card_max);
    fmpz_t m;fmpz_init_set_ui(m,1);
    for(;fmpz_cmp_ui(m,card_max)<=0;fmpz_add_ui(m,m,1)){
        fmpz_gcd(A,m,n);
        if(fmpz_cmp_ui(A,1)>0){
            fmpz_clear(a);fmpz_clear(A);fmpz_clear(m);
            return false;
        }
    }
    /******************************Step 4*******************************/
    if(t>k+1)
        card_max=t-k-1;
    else
        card_max=0;
    card_max=exp2l(card_max);
    /** Declaration of polys in Z/nZ **/
    fmpz_mod_poly_t poly1,poly2,modulo,poly2_core,polytmp,poly_one;
    /** Set poly1 **/
    fmpz_mod_poly_init(poly1,n);
    fmpz_mod_poly_set_coeff_ui(poly1,0,1);
    /** Set poly2 **/
    fmpz_mod_poly_init(poly2,n);

    /** Set modulo **/
    fmpz_mod_poly_init(modulo,n);
    fmpz_mod_poly_set_coeff_fmpz(modulo,0,a);
    ulong exp2_t=exp2l(t);
    fmpz_mod_poly_set_coeff_ui(modulo,2*exp2_t,1);
    fmpz_sub_ui(A,n,2);
    fmpz_mod_poly_set_coeff_fmpz(modulo,exp2_t,A);
    /** Set poly2_core **/
    fmpz_mod_poly_init(poly2_core,n);
    fmpz_mod_poly_set_coeff_ui(poly2_core,exp2_t,2);
    fmpz_mul_si(A,a,-1);
    fmpz_mod_poly_set_coeff_fmpz(poly2_core,0,A);
    fmpz_fdiv_q_2exp(A,n,t+1);
    fmpz_mod_poly_powmod_fmpz_binexp(poly2_core,poly2_core,A,modulo);
    fmpz_fdiv_r_2exp(A,n,t+1);temp_ui=fmpz_get_ui(A);
    fmpz_mod_poly_shift_left(poly2_core,poly2_core,temp_ui);
    /** Set polytmp **/
    fmpz_mod_poly_init(polytmp,n);
    /** Set poly_one **/
    fmpz_mod_poly_init(poly_one,n);
    fmpz_mod_poly_set_coeff_ui(poly_one,0,1);
    fmpz_set_ui(m,1);
    for(;fmpz_cmp_ui(m,card_max)<=0;fmpz_add_ui(m,m,1)){
        fmpz_mod_poly_set_coeff_fmpz(poly1,1,m);
        fmpz_mod_poly_powmod_fmpz_binexp(polytmp,poly1,n,modulo);
        fmpz_mod_poly_scalar_mul_fmpz(poly2,poly2_core,m);
        fmpz_mod_poly_add(poly2,poly2,poly_one);
        if(fmpz_mod_poly_equal(polytmp,poly2)){
            fmpz_clear(a);
            fmpz_clear(A);
            fmpz_clear(m);
            fmpz_mod_poly_clear(poly1);
            fmpz_mod_poly_clear(poly2);
            fmpz_mod_poly_clear(poly2_core);
            fmpz_mod_poly_clear(modulo);
            fmpz_mod_poly_clear(polytmp);
            fmpz_mod_poly_clear(poly_one);
            return false;
        }
    }
    fmpz_clear(a);
    fmpz_clear(A);
    fmpz_clear(m);
    fmpz_mod_poly_clear(poly1);
    fmpz_mod_poly_clear(poly2);
    fmpz_mod_poly_clear(poly2_core);
    fmpz_mod_poly_clear(modulo);
    fmpz_mod_poly_clear(polytmp);
    fmpz_mod_poly_clear(poly_one);
    return true;
}
bool AKS_Berrizbeitia(fmpz_t n){
    if(fmpz_cmp_ui(n,101)<0)//Check the small numbers
        return PRIME_TABLE[fmpz_get_ui(n)];
    fmpz_t small_p_mul;fmpz_init(small_p_mul);
    fmpz_set_str(small_p_mul,"2305567963945518424753102147331756070",10);
    fmpz_gcd(small_p_mul,small_p_mul,n);//Trial division for small prime numbers using Euclid method
    if(fmpz_cmp_ui(small_p_mul,1)>1){
        fmpz_clear(small_p_mul);
        return false;
    }
    mpz_t mpz_n;mpz_init(mpz_n);
    fmpz_get_mpz(mpz_n,n);
    if(mpz_perfect_power_p(mpz_n)){//Check the perfect power
        mpz_clear(mpz_n);
        return false;
    }
    mpz_clear(mpz_n);
    fmpz_mod_ui(small_p_mul,n,4);
    if(fmpz_get_ui(small_p_mul)==1){//The case of n=4k+1
        fmpz_clear(small_p_mul);
        return CASE_I(n);
    }
    fmpz_clear(small_p_mul);
    return CASE_II(n);//The case of n=4k+3
}
#include<iostream>
using namespace std;
int main(){
    long long unsigned k=exp2l(64)-353;
    cout<<k<<endl;
    fmpz_t e;fmpz_init_set_ui(e,k);
    cout<<AKS_Berrizbeitia(e)<<endl;
    return 0;
}
