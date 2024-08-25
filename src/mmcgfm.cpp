// This script implement multi-study multi-modality covriate-augmented generalized factor model based on variational inference.
// Date: 2024-06-06

// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
 //' @noRd
 //'
 // diag(W0* Cki * W0^t)
 vec decomp(const mat& Cki, const mat& W0){
   vec s, tmp1;
   mat U, V, WC12;
   svd(U, s, V, Cki);
   WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
   tmp1 = sum(WC12 % WC12, 1);
   return tmp1;
 }

void update_betaf(const field<mat>& Xf, const field<vec>& tauf, const field<mat>& Zf, const field<mat>& Af, const field<mat>& Bf,
                  const field<mat>&  Mf, const field<mat>&  Of, const field<vec>&  wf, field<mat>& betaf){
  int s, id, im, pm, dz = Zf(0).n_cols, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  mat zzt(dz, dz, fill::zeros);
  for(s=0; s<S; ++s){ // Loop for all sources
      zzt += Zf(s).t()*Zf(s);
  }
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
      pm = Xf(im,id, 0).n_cols;
      mat mat_tmp = zeros(dz, pm);
      for(s=0; s<S; ++s){ // Loop for all sources
        if(Bf(im,id, s).n_cols>1){
          mat_tmp += Zf(s).t()*(Xf(im,id, s) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm) -  Mf(s)*Af(im,id).t() - Of(s)*Bf(im, id, s).t()); //dz*p_m
        }
      }
      betaf(im, id) = zzt.i()* mat_tmp; // dz*p_m
    }
  }
  
}

void update_Af_fast(const field<mat>& Xf, const field<vec>& tauf, const field<mat>& Zf, const field<mat>& Bf,const field<mat>&  Mf, 
               const cube& Sigma, const field<mat>&  Of, const field<vec>&  wf, const field<mat>& betaf,
               field<mat>& Af){ // Whether weighted using invLambdaf? by performance
  
  int s, j, ns, id, im, pm, q = Mf(0).n_cols, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  mat mat_tmp, mat_tmpx;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
      pm = Xf(im,id, 0).n_cols;
      mat_tmp = zeros(q, q);
      mat_tmpx = zeros(pm, q);
      if(Bf(im,id, 0).n_cols>1){ // If the modality exists, do .
        for(s=0; s<S; ++s){ // Loop for all sources
          ns = Xf(0,0,s).n_rows;
          mat_tmp += Mf(s).t()*Mf(s) + ns*Sigma.slice(s); // q*q
          mat mat_tmpL = Xf(im, id, s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  - Of(s)*Bf(im, id, s).t(); // ns * p_m
          mat_tmpx += trans(mat_tmpL) * Mf(s); //p_m*q
        }
        Af(im, id) = mat_tmpx * mat_tmp.i(); // pm*q
      }
    }
  }
}


void update_Af(const field<mat>& Xf, const field<vec>& tauf, const field<mat>& Zf, const field<mat>& Bf,const field<mat>&  Mf, 
               const cube& Sigma, const field<mat>&  Of, const field<vec>&  wf, const field<mat>& betaf,
               const field<vec>&  invLambdaf, field<mat>& Af){ // Whether weighted using invLambdaf? by performance
  
  int s, j, ns, id, im, pm, q = Mf(0).n_cols, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  mat mat_tmp, mat_tmpx;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
      pm = Xf(im,id, 0).n_cols;
      mat_tmpx = zeros(pm, q);
      if(Bf(im,id, 0).n_cols>1){
        for(s=0; s<S; ++s){ // Loop for all sources
         ns = Xf(0,0,s).n_rows;
          // mat_tmp += Mf(s).t()*Mf(s) + ns*Sigma.slice(s); // q*q
          mat mat_tmpL = Xf(im, id, s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  - Of(s)*Bf(im, id, s).t(); // ns * p_m
          mat_tmpx += trans(mat_tmpL % repmat(invLambdaf(im, id, s).t(), ns, 1)) * Mf(s); //p_m*q
        }
        for(j =0; j<pm; ++j){
          mat_tmp = zeros(q, q);
          for(s=0; s<S; ++s){ // Loop for all sources
            ns = Xf(0,0,s).n_rows;
            mat_tmp += (Mf(s).t()*Mf(s) + ns*Sigma.slice(s))*invLambdaf(im, id, s)(j); // q*q
            
          }
         Af(im, id).row(j) = mat_tmpx.row(j) * mat_tmp.i(); // pm*q
          
        }
      }
    }
  }
}




void update_Bf(const field<mat>& Xf, const field<vec>& tauf, const field<mat>& Zf, const field<mat>& Af,const field<mat>&  Mf, 
      const field<mat>& Phi, const field<mat>&  Of, const field<vec>&  wf, const field<mat>& betaf, field<mat>& Bf){
  int s, ns, id, im, pm, qs, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  mat mat_tmp;
  for(s=0; s<S; ++s){ // Loop for all sources
   qs = Of(s).n_cols;
   ns = Xf(0,0,s).n_rows;
   for(id=0; id <d; ++id){ // Loop for all modality group
     for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
      pm = Xf(im,id, 0).n_cols;
      if(Bf(im,id, s).n_cols>1){
          mat_tmp = Of(s).t()*Of(s) + ns*Phi(s); // qs*qs
          mat mat_tmpx = trans(Xf(im, id, s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  - Mf(s)*Af(im, id).t()) * Of(s); //p_m*qs
          Bf(im, id, s) = mat_tmpx * mat_tmp.i(); // pm*qs
       }
      }
      
    }
  }
}

void update_invLambdaf(const field<mat>& Xf, const field<vec>& tauf, const field<mat>& Zf,const field<mat>& betaf,  const field<mat>& Af, 
                       const field<mat>& Bf, const field<mat> Sf_y, 
                       const field<mat>&  Mf, const cube& Sigma, const field<mat>&  Of, const field<mat>&  Phi,const field<vec>&  wf, 
                       const cube& zeta, field<vec>&  invLambdaf){
  
  int s, id, im, pm, ns,  d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  for(s=0; s<S; ++s){ // Loop for all sources
    ns = Xf(0,0, s).n_rows;
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        pm = Xf(im,id, 0).n_cols;
        if(Bf(im,id, s).n_cols>1){
          vec tmp_vec1 = decomp(Sigma.slice(s), Af(im, id)); // 
          vec tmp_vec2 = decomp(Phi(s), Bf(im, id, s)); // 
          mat dX = (Xf(im, id, s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  - Mf(s)*Af(im, id).t() -Of(s)*Bf(im, id, s).t()) ; //p_m*qs
          tmp_vec1 = (tmp_vec1 + tmp_vec2) + trans(mean(dX % dX + Sf_y(im,id,s))) +  zeta(im, id, s);
          invLambdaf(im, id, s) = 1.0/tmp_vec1; // pm*1
        }
      }
      
    }
  }
}




void update_sigma2(const field<vec>&  wf, const cube& zeta, const field<mat>& Bf,  cube& sigma2){
  int s, id, im, pm, ns,  d = wf.n_cols, c=wf.n_rows, S = wf.n_slices;
  mat mat_tmp;
  for(s=0; s<S; ++s){ // Loop for all sources
    ns = wf(0,0, s).n_elem;
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id, s).n_cols>1){
          
          sigma2(im, id, s) = accu(wf(im, id, s) % wf(im, id, s))/ns +zeta(im, id, s);
        }
      }
      
    }
  }
}


void add_IC_Orth(field<mat>& Af, field<mat>& Bf){
  // Try the orthogonal matrix method
  int id, im, s, qs1, q=Af(0,0).n_cols, d = Bf.n_cols, c=Bf.n_rows, S=Bf.n_slices;
  
  
  for(s=0; s<S; ++s){
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id, s).n_cols>1){
          qs1 = Bf(im,id, s).n_cols;
          if(s==1){
            mat AB1 = join_rows(Af(im, id), Bf(im, id, s)); // pm * (q+q1)
            mat U1, V1;
            vec s1;
            svd(U1, s1, V1, AB1);
            vec signU1 = sign(U1.row(0).t());
            Af(im, id) = U1.cols(0,q-1) * diagmat(s1.subvec(0, q-1) % signU1.subvec(0,q-1));
            Bf(im, id, s) = U1.cols(q, q+qs1-1) * diagmat(s1.subvec(q, q+qs1-1) % signU1.subvec(q, q+qs1-1));
          }else{
            mat U1, V1;
            vec s1;
            svd(U1, s1, V1, Bf(im, id, s));
            vec signU1 = sign(U1.row(0).t());
            Bf(im, id, s) = U1.cols(0,qs1-1) * diagmat(s1.subvec(0, qs1-1) % signU1.subvec(0,qs1-1));
          }
          
        }
      }
    }
  }
}


// void add_IC_beta(const mat&Z, const field<mat>& Bf, field<vec>& betaf, field<vec>& Xif){
//   
//   int id, im, qs1,  d = Bf.n_cols, c=Bf.n_rows;
//   vec vec_tmp;
//   for(id=0; id <d; ++id){ // Loop for all modality group
//     for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
//       if(Bf(im,id).n_cols>1){
//         vec_tmp = (Z.t()*Z).i() * Z.t() * Xif(im, id);
//         betaf(im,id) += vec_tmp;
//         Xif(im, id) -= Z*vec_tmp;
//       }
//     }
//   }
//   
// }
double ELBOfun(const field<mat>& Xf, const vec& typeID, const field<vec>& tauf, const field<mat>& Zf,const field<mat>& betaf,
               const field<mat>& Af, const field<mat>& Bf, const  field<vec>&  invLambdaf, const cube& sigma2,
               const field<mat>& Muf_y, const field<mat>& Sf_y,
               const field<mat>&  Mf, const cube& Sigma, const field<mat>&  Of, const field<mat>& Phi, const field<vec>&  wf,
               const cube&  zeta){
  // Basic information:
  int s, ns, id, im, pm, qs, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  double tmp_v;
  
  // log P(X|Y) 
  double val_tmp0 = 0.0;
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    for(id=0; id <d; ++id){ // Loop for all modality group
      if(typeID(id)==2){
        for(im=0; im<c; ++im){ // loop for all modalities in Count-type group.
          if(Bf(im,id,s).n_cols>1){ 
            val_tmp0 += accu(Xf(im,id,s) % Muf_y(im,id,s)-exp(Muf_y(im,id, s) + Sf_y(im,id, s)/2));
          }
        }
      }
      if(typeID(id)==3){
        for(im=0; im<c; ++im){ // loop for all modalities in Binomial-type group.
          if(Bf(im,id, s).n_cols>1){ 
            rowvec nvec = max(Xf(im,id, s));
            val_tmp0 += accu((Xf(im,id, s)-repmat(nvec, ns, 1)) % Muf_y(im,id, s));  //
            val_tmp0 += -accu(repmat(nvec, ns, 1) %  exp(-Muf_y(im,id, s) + Sf_y(im,id,s)/2)); // term2 is not written now!
          }
        }
      }
    }
    
  }
  
  
  // log P(Y|F, H, Z)
  //Rprintf("Good EBLO1!\n");
  double val_tmp1 = 0.0; //
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id, s).n_cols>1){ 
          pm = Xf(im,id, s).n_cols;
          mat mat_tmp = (Muf_y(im,id,s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  - Mf(s)*Af(im, id).t() - Of(s)*Bf(im, id, s).t());
          vec vec_tmp = ns*decomp(Sigma.slice(s), Af(im,id)) + ns*decomp(Phi(s), Bf(im,id,s))  + trans(sum(Sf_y(im,id,s))) + ns* zeta(im,id, s);
          val_tmp1 += accu(vec_tmp % invLambdaf(im,id, s)) + accu(sum(mat_tmp % mat_tmp) % invLambdaf(im,id, s).t()) - ns*accu(log(invLambdaf(im,id, s)));
        }
      }
    }
  }
  val_tmp1 = -0.5* val_tmp1;
  
  //Rprintf("Good EBLO2!\n");
  // log P(F ) + log P(H) + log P(V)
  double val_tmp2 = 0.0, val_tmp3 = 0.0;
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id,s).n_cols>1){ 
          tmp_v = log(sigma2(im,id,s)+1e-10);
          val_tmp2 += ns * tmp_v + accu(wf(im,id, s) % wf(im,id, s) + zeta(im,id, s))/sigma2(im,id,s);
        }
      }
    }
    val_tmp3 +=  accu(Mf(s) % Mf(s)) + ns*trace(Sigma.slice(s)) + accu(Of(s) % Of(s)) + ns*trace(Phi(s)) ;
  }
  val_tmp2 = -0.5* val_tmp2;
  val_tmp3 = -0.5* val_tmp3;
  
  
  // Entropy
  double entropy= 0.0;
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    entropy += ns* log(det(Sigma.slice(s))) + ns* log(det(Phi(s)));
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id).n_cols>1){ 
          tmp_v = log(zeta(im,id, s));
          entropy += ns*tmp_v + accu(log(Sf_y(im,id, s)+1e-10));
          
        }
      }
    }
  }
  entropy = 0.5*entropy;
  //Rprintf("Good EBLO4!\n");
  
  return val_tmp0+val_tmp1+val_tmp2+val_tmp3+ entropy;
}


void VB_Estep(const field<mat>& Xf, const vec& typeID, const field<vec>& tauf, const field<mat>& Zf,const field<mat>& betaf,
              const field<mat>& Af, const field<mat>& Bf, const  field<vec>&  invLambdaf, const cube& sigma2,
              field<mat>& Muf_y, field<mat>& Sf_y,
               field<mat>&  Mf, cube& Sigma, field<mat>&  Of, field<mat>& Phi, field<vec>&  wf, cube&  zeta){
  
  // Basic information:
  int s, ns, id, im, pm, qs, q = Af(0,0).n_cols, d = Xf.n_cols, c=Xf.n_rows, S = Xf.n_slices;
  // Rprintf("Good Estep1!\n");
  // double elbo0 =  ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
  // Update Muf_y and Sf_y
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    for(id=0; id <d; ++id){ // Loop for all modality group
      // Rprintf("Good Estep: s = %d, id = %d!\n", s, id);
      // Initialize in the external if typeID(id)==1
      // Rprintf("Good Estep1.1!\n");
      if(typeID(id)==2){
        for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
          if(Bf(im,id,s).n_cols>1){
            pm = Xf(im,id,s).n_cols;
            mat Mu_y1 = Muf_y(im,id, s);
            mat tZ = Zf(s)*betaf(im, id) + repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  + Mf(s)*Af(im, id).t() + Of(s)*Bf(im, id, s).t();
            Muf_y(im,id, s) = (Xf(im,id, s) -  exp(Mu_y1) % (1-Mu_y1) + repmat(invLambdaf(im,id, s).t(), ns, 1)% tZ) /
              (exp(Mu_y1) +  repmat(invLambdaf(im,id, s).t(), ns, 1) );
            Sf_y(im,id, s)= 1.0 / (exp(Muf_y(im,id,s)) + repmat(invLambdaf(im,id,s).t(), ns, 1) );
          }
        }
      }
      if(typeID(id)==3){
        for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
          if(Bf(im,id, s).n_cols>1){
            pm = Xf(im,id, s).n_cols;
            rowvec nvec = max(Xf(im,id, s));
            mat tZ = Zf(s)*betaf(im, id) + repmat(tauf(im,id, s)+wf(im, id, s),1, pm)  + Mf(s)*Af(im, id).t() + Of(s)*Bf(im, id, s).t();
            mat Mu_y2 = Muf_y(im,id, s); // take out this submatrix
            Mu_y2 = (Xf(im,id, s)- (1/(1+exp(-Mu_y2))) % repmat(nvec,ns, 1) + repmat(invLambdaf(im,id, s).t(), ns, 1)% tZ) /
              ((1/(1+exp(-Mu_y2))) % repmat(nvec,ns, 1) + repmat(invLambdaf(im,id, s).t(), ns, 1) );
            Muf_y(im,id, s) = Mu_y2;
            Sf_y(im,id, s) = 1.0 / ( (1/(1+exp(-Mu_y2))) %(1- 1/(1+exp(-Mu_y2))) % repmat(nvec,ns, 1) + repmat(invLambdaf(im,id, s).t(), ns, 1));
          }
        }
      }
    }
  }
  
  // double elbo1 =  ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
  // Rprintf("dMu_y= %4f\n", elbo1 - elbo0);
  // Update Si and mi
  mat Sigma_inv(q,q, fill::zeros);
  double val_tmp1;
  // Rprintf("Good Estep1!\n");
  //mat mat_tmp3;
  //field<mat> Om_inv(c,d);
  for(s=0; s<S; ++s){
    ns = Xf(0,0, s).n_rows;
    qs = Bf(0,0,s).n_cols;
    mat Phi_inv(qs, qs, fill::zeros);
    Sigma_inv = zeros(q,q);
    mat M_tmp(ns,q, fill::zeros);
    mat O_tmp(ns,qs, fill::zeros);
    double elbo11 =  ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<c; ++im){ // loop for all modalities in each modality group.
        // Rprintf("Good Estep: s = %d, id = %d, im = %d!\n", s, id, im);
        if(Bf(im,id, s).n_cols>1){
          pm = Xf(im,id, s).n_cols;
          mat AfL = (Af(im,id) % repmat(invLambdaf(im,id, s), 1, q));
          mat BfL = (Bf(im,id,s) % repmat(invLambdaf(im,id, s), 1, qs));
          Sigma_inv += Af(im,id).t() * AfL;
          Phi_inv += Bf(im,id,s).t() * BfL;
          // Rprintf("Good Estep2: s = %d, id = %d, im = %d!\n", s, id, im);
          
          val_tmp1 = accu(invLambdaf(im,id, s));
          // M_tmp += (Xf(im,id) - repmat(Af(im,id),1, pm) - repmat(muf(im,id).t(), n, 1)-Xif(im,id)*Bf(im,id).t())*BfL;
          mat mat_tmp2 = Muf_y(im,id,s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s),1, pm)  - Mf(s)*Af(im, id).t() - Of(s)*Bf(im, id, s).t();
          //Om_inv(im,id) = mat_tmp1 + diagmat(1.0/Sigmamf(im,id));
          //Of(im,id) = Om_inv(im,id).i();
          zeta(im,id, s) = 1/(val_tmp1 + 1.0/sigma2(im,id, s)); // zetaf is not related to i
          // mat_tmp3 = Om_inv(im,id).i(); // Only update Xif
          //Xif(im,id) = (mat_tmp2 - (M*Bf(im,id).t())*BfL)* mat_tmp3;
          wf(im,id, s) =zeta(im,id, s) * sum(mat_tmp2 % repmat(invLambdaf(im,id, s).t(),ns , 1), 1); // sum by row
          
          // Rprintf("Good Estep3: s = %d, id = %d, im = %d!\n", s, id, im);
          mat mat_tmp1 = (Muf_y(im,id,s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im,id,s),1, pm) -Of(s)*Bf(im, id, s).t())*AfL;
          mat mat_tmp3 = (Muf_y(im,id,s) - Zf(s)*betaf(im, id) - repmat(tauf(im,id, s)+wf(im,id,s),1, pm) -Mf(s)*Af(im, id).t())*BfL;
          M_tmp += mat_tmp1; // Must first update Xif, then update M_tmp
          O_tmp += mat_tmp3;
          
        }
      }
      
    }
    // double elbo21 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dMO_w_zeta= %4f, s=%d \n", elbo21 - elbo11, s);
    // Rprintf("Out S: Good Estep4: s = %d!\n", s);
    Sigma.slice(s) = (Sigma_inv + eye(q,q)).i(); // this is some problem
    // Phi_inv.print();
    // Phi.slice(s).print();
    Mf(s) = M_tmp * Sigma.slice(s);
    // double elbo2 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dM_Sigma= %4f, s=%d \n", elbo2 - elbo21, s);
    
    Phi(s) = (Phi_inv+eye(qs, qs)).i();
    Of(s) = O_tmp * Phi(s);
    // double elbo22 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dO_Phi= %4f, s=%d \n", elbo22 - elbo2, s); 
  }
  
  
  
}


// [[Rcpp::export]]
Rcpp::List vb_mmgfmcpp(const Rcpp::List& XList, const arma::vec& typeID,
                       const arma::mat& numvarmat, const Rcpp::List&  tauList,
                       const Rcpp::List& Zlist, const Rcpp::List& betalist_int,
                       const Rcpp::List& Alist_int, const Rcpp::List& Blist_int,
                       const Rcpp::List& invLambdalist_int, const arma::cube& sigma2_int,
                       const Rcpp::List& Mulist_y_int, const Rcpp::List& Slist_y_int, // improve by puting them into cpp not in R.
                       const Rcpp::List& Mlist_int, const arma::cube& Sigma_int,
                       const Rcpp::List& Olist_int, const Rcpp::List& Philist_int, 
                       const arma::cube& zeta_int, const Rcpp::List&  wlist_int, 
                       const double& epsELBO, const int& maxIter, const bool& verbose, const bool& A_fast,
                       const bool& add_IC_iter=true){
  // XList is a embeded list, in the first level: the data from S sources; in the second level: the modalities with different types
  // Alist_int is a list, each entry is a p_m * q matrix; Blist_int is a embeded list, in second level, a p_m* q_s matrix.
  // typeID: 1 means Gaussian; 2 means Poisson; 3 means Binomial;
  // Xf field represents the modality group with same variable type.
  // numvarmat: the number ofvariables in each modality within modality group.
  // betalist_int[[t]] \in R^(d*M_t)
  
  int i, im, s, S = XList.length(), d = numvarmat.n_rows,  c= numvarmat.n_cols; 
  // where c*d = M, the total number of modalities.
  // S is the number of sources; d is the number of variable types; c is the max number of modalities with same type.
  int d2 = typeID.n_elem;
  if(d != d2){
    stop("The length of XList must be equal to the length of typeID!");
  }
  // Rprintf("Good entry1!\n");
  // Rprintf("d=%d, c=%d, S=%d\n", d, c, S);
  
//  XList: Xf
//  tauList: tauf
//  betalist_int: betaf
//  Alist_int: Af
//  Blist_int: Bf
// invLambdalist_int: invLambdaf
// sigma2_int:  sigma2
//   Mulist_y_int: Muf_y
//   Slist_y_int: Sf_y
//   Mlist_int: Mf
//   Sigma_int: Sigma
//  Olist_int: Of
// Phi_int: Phi
// zetalist_int: zetaf
// wlist_int: wf
  
  
  field<mat> Xf(c,d,S),  Muf_y(c,d,S),  Sf_y(c,d,S), Bf(c,d,S),  Mf(S), Of(S), Zf(S), Phi(S); // 
  // Xf, Muf_y and Sf_y's entry is a n_s * p_(c,d) matrix; Bf's entry is p_(c,d)*q_s
  // tauf is the c*d filed<cube>, each entry is a n_s vector, the offset.
  field<vec> tauf(c,d, S), wf(c,d,S), invLambdaf(c,d, S);


  
  cube sigma2(sigma2_int), zeta(zeta_int);
  cube Sigma(Sigma_int);
  Rprintf("Initialize the paramters related to s!\n");
  // Initialize the paramters related to s
  for(s=0; s<S; ++s){
    // Rprintf("s= %d\n", s);
    mat Mtmp = Mlist_int[s];
    Mf(s) = Mtmp;
    mat Otmp = Olist_int[s];
    Of(s) = Otmp;
    mat Ztmp = Zlist[s];
    Zf(s) = Ztmp;
    mat phitmp = Philist_int[s];
    Phi(s) = phitmp;
    Rcpp::List Xtmplist = XList[s];
    Rcpp::List Muytmplist = Mulist_y_int[s];
    Rcpp::List Sytmplist = Slist_y_int[s];
    Rcpp::List Btmplist = Blist_int[s];
    Rcpp::List tmptaulist = tauList[s];
    Rcpp::List tmpwlist = wlist_int[s];
    Rcpp::List tmpinvLamlist = invLambdalist_int[s];
    // Rprintf("Good entry2!\n");
    for(i=0; i<d; ++i){ // put the modality group matrix of list into a field.
      // Rprintf("s= %d, i = %d \n", s, i);
      mat Xtmp = Xtmplist[i];
      mat Mutmp = Muytmplist[i];
      mat Stmp = Sytmplist[i];
      // Rprintf("s= %d, i = %d \n", s, i);
      mat Btmp = Btmplist[i];
      mat tautmp = tmptaulist[i];
      vec invLamtmp = tmpinvLamlist[i];
      // Rprintf("s= %d, i = %d \n", s, i);
      mat wtmp = tmpwlist[i];
      // Rprintf("Good entry3!\n");
      vec p_vec(c+1, fill::zeros); // initilize
      for(im=0; im<c; ++im){
        if(numvarmat(i,im)>0){
          // Rprintf("s= %d, i = %d, im=%d \n", s, i, im);
          
          p_vec(im+1) = p_vec(im) + numvarmat(i,im);
          Xf(im,i, s) = Xtmp.cols(p_vec(im), p_vec(im+1)-1);
          Muf_y(im,i, s) = Mutmp.cols(p_vec(im), p_vec(im+1)-1);
          Sf_y(im,i, s) = Stmp.cols(p_vec(im), p_vec(im+1)-1);
          
          // Rprintf("Good entry1!\n");
          Bf(im,i, s) = Btmp.rows(p_vec(im), p_vec(im+1)-1);
          tauf(im,i, s) = tautmp.col(im);
          invLambdaf(im,i, s) = invLamtmp.subvec(p_vec(im), p_vec(im+1)-1);
          wf(im, i, s) = wtmp.col(im); // is a n_s vector.
        }else{
          Xf(im,i, s) = zeros(1,1); // fill 1-by-1 zero matrix for other empty position.
          Muf_y(im,i, s) = Xf(im,i,s);
          Sf_y(im,i, s) = Xf(im,i, s);
          Bf(im,i,s) = Xf(im,i, s);
          tauf(im,i, s) =zeros(1,1);
          invLambdaf(im,i, s) = Xf(im,i, s);
          wf(im, i, s) = Xf(im,i, s);
        }
        
      }
    }
    
  }
  
  Rprintf("Initialize the paramters not related to s!\n");
  field<mat> Af(c,d), betaf(c,d);
  // Initialize the paramters not related to s
  for(i=0; i<d; ++i){ // put the modality group matrix of list into a field.
    
    mat Atmp = Alist_int[i];
    mat bbtmp = betalist_int[i];
    vec p_vec(c+1, fill::zeros); // initilize
    for(im=0; im<c; ++im){
      if(numvarmat(i,im)>0){
        // Rprintf("s= %d, i = %d, im=%d \n", s, i, im);
        
        p_vec(im+1) = p_vec(im) + numvarmat(i,im);
        betaf(im,i) = bbtmp.cols(p_vec(im), p_vec(im+1)-1);
        // Rprintf("Good entry1.1!\n");
        Af(im,i) = Atmp.rows(p_vec(im), p_vec(im+1)-1);
      }else{
        betaf(im,i) = zeros(1,1);
        Af(im,i) = zeros(1,1);
      }
      
    }
  }
  
  
  
  
  // Rprintf("Good entry2!\n");
  // Initialize
  
  vec ELBO_vec(maxIter);
  ELBO_vec(0) = -1e20;
  int iter;
  
  for(iter = 1; iter < maxIter; ++iter){
    
    
    // Rprintf("E step starting!\n");
    // double elbo0 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("elbo0=%4f!\n", elbo0);
    // VB E-step
    VB_Estep( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y,  Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("Finish E step!\n");
    //VB M-step
    
    // double elbo1 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    
    
    // update mu
    // Rprintf("update beta\n");
    update_betaf(Muf_y, tauf, Zf, Af, Bf, Mf, Of, wf, betaf);
    // double elbo2 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dbeta= %4f \n", elbo2 - elbo1);
    
    // update beta
    // Rprintf("update A\n"); // This is some problem
    if(A_fast){
      update_Af_fast( Muf_y, tauf, Zf, Bf, Mf, Sigma, Of, wf, betaf, Af);
    }else{
      update_Af( Muf_y, tauf, Zf, Bf, Mf, Sigma, Of, wf, betaf,invLambdaf, Af);
    }
    
    // double elbo21 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dA= %4f \n", elbo21 - elbo2);
    
    //update B
    // Rprintf("update B\n");//there is some problem in updating lambda
    update_Bf(Muf_y, tauf, Zf, Af, Mf, Phi, Of, wf, betaf, Bf) ;
    if(add_IC_iter){
      add_IC_Orth(Af, Bf); // Add identifiability condition
    }
    
    // double elbo3 = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dB= %4f \n", elbo3 - elbo21);
    
    // update Lambda
    // Rprintf("update Lambda\n"); // 
    update_invLambdaf( Muf_y, tauf, Zf, betaf, Af, Bf, Sf_y, Mf, Sigma, Of, Phi, wf, zeta, invLambdaf);
    // double elbo4 =   ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dLambda= %4f\n", elbo4 - elbo3);
    
    // update Sigmam
    
    // Rprintf("update Sigmam\n");
    update_sigma2( wf, zeta, Bf, sigma2);
    // double elbo5 =  ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    // Rprintf("dsigma= %4f\n", elbo5 - elbo4);
    
    
    
    
    ELBO_vec(iter) = ELBOfun( Xf, typeID, tauf, Zf, betaf, Af, Bf, invLambdaf, sigma2, Muf_y, Sf_y, Mf, Sigma, Of, Phi, wf, zeta);
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n",
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  if(!add_IC_iter){
    add_IC_Orth(Af, Bf); // Add identifiability condition
    // VB_Estep(Xf, Af, typeID, Z, Muf_y, Sf_y, M, S, Xif,  Om, Bf, betaf, Sigmam, muf, invLambdaf);
  }
  //add_IC_beta(Z, Bf, betaf, Xif);
  // Put Bf into a list
  Rcpp::List BList(S), invLambdaList(S), WList(S);
  field<mat> tmpB;
  field<vec> tmpLam, tmpw;
  for(s=0; s<S; ++s){
    tmpB = Bf.slice(s);
    BList[s] = tmpB;
    tmpLam = invLambdaf.slice(s);
    invLambdaList[s] = tmpLam;
    tmpw = wf.slice(s);
     WList[s] = tmpw;
  }
  
  
  // output return value
  List resList = List::create(
    Rcpp::Named("hbeta") = betaf,
    Rcpp::Named("hA") = Af,
    Rcpp::Named("hB") = BList, // Put the Bf into a list?????, then ouput????
    Rcpp::Named("hF") = Mf,
    Rcpp::Named("hH") = Of,
    Rcpp::Named("hSigma") = Sigma,
    Rcpp::Named("hPhi") = Phi,
    Rcpp::Named("hv") = WList,
    Rcpp::Named("hzeta") = zeta,
    Rcpp::Named("hsigma2") = sigma2, // output is not
    Rcpp::Named("hinvLambda") = invLambdaList,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}

