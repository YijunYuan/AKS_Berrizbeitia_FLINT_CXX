#ifndef IQF_CALC_H_INCLUDED
#define IQF_CALC_H_INCLUDED
#include<flint/fmpz.h>
void inline IQF_MUL    (           fmpz_t base,fmpz_t n,fmpz_t u1,fmpz_t v1,fmpz_t u2,fmpz_t v2,fmpz_t uo,fmpz_t vo){
    fmpz_t temp;fmpz_init_set(temp,base);
    fmpz_mul(temp,temp,v1  );
    fmpz_mul(temp,temp,v2  );
    fmpz_set(uo  ,u1       );
    fmpz_mul(uo  ,uo  ,u2  );
    fmpz_add(uo  ,uo  ,temp);
    fmpz_mod(uo  ,uo  ,n   );
    fmpz_set(temp,u1       );//temp=u1
    fmpz_mul(temp,temp,v2  );//temp=u1*v2
    fmpz_set(vo  ,u2       );//vo=u2
    fmpz_mul(vo  ,vo  ,v1  );//vo=u2*v1
    fmpz_add(vo  ,vo  ,temp);
    fmpz_mod(vo  ,vo  ,n   );
    fmpz_clear(temp);
    return ;
};//(uo,vo)=(u1,v1)*(u2,v2)
void inline IQF_SQRT   (           fmpz_t base,fmpz_t n,fmpz_t u1,fmpz_t v1,                    fmpz_t uo,fmpz_t vo){
    IQF_MUL(base,n,u1,v1,u1,v1,uo,vo);
    return ;
};//(uo,vo)=(u1,v1)^2
void inline IQF_EXP    (fmpz_t exp,fmpz_t base,fmpz_t n,fmpz_t u1,fmpz_t v1,                    fmpz_t uo,fmpz_t vo){
    fmpz_t t;fmpz_init_set(t,exp);
    fmpz_set_ui(uo,1);fmpz_set_ui(vo,0);
    fmpz_t w1,w2;fmpz_init_set(w1,u1);fmpz_set(w2,v1);
    fmpz_t t1,t2;fmpz_init(t1);fmpz_init(t2);
    while(fmpz_is_zero(t)==0){
        if(fmpz_is_odd(t)){
            fmpz_set(t1,uo);fmpz_set(t2,vo);
            IQF_MUL(base,n,t1,t2,w1,w2,uo,vo);
        }
        fmpz_fdiv_q_2exp(t,t,1);
        fmpz_set(t1,w1);fmpz_set(t2,w2);
        IQF_SQRT(base,n,t1,t2,w1,w2);
    }
    fmpz_clear(t);fmpz_clear(w1);fmpz_clear(w2);fmpz_clear(t1);fmpz_clear(t2);
    return ;
};//(uo,vo)=(u1,v1)^p
#endif // IQF_CALC_H_INCLUDED
