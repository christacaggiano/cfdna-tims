library("data.table")
library("bigstatsr")
library("dplyr")
library("pROC")

########################################################
# separate data into input, outcome, and covariates 

meths = fread("../data/final_case_control_both_sites.txt", data.table=F)
status = fread("../data/final_case_control_both_status.txt", data.table=F)

covs = meths[, 1:8]
meth_only = meths[, -c(1:8)]

ucsf_num = nrow(status[status$group == 0, ])

ucsf = meth_only[1:ucsf_num, ]

ucsf_covs = covs[1:ucsf_num, c("sex female", "age", "concentration", "starting", "race white")]
ucsf_status = status[1:ucsf_num, c("sample_type") ]

uq = meth_only[(ucsf_num+1):nrow(meth_only), ]
uq_status = status[(ucsf_num+1):nrow(meth_only), c("sample_type")]

uq_covs = covs[(ucsf_num+1):nrow(meth_only),c("sex female", "age", "concentration", "starting", "race white")]

print(dim(ucsf))
print(dim(uq))

########################################################
##### UCSF trained UQ tested

system("rm ~/temp_matrix.bk")
fitdat=as_FBM(scale(ucsf), backingfile = "~/temp_matrix")

res.all.ucsf=big_spLogReg(X=fitdat, y01.train=ucsf_status, alphas=1e-4, , 
                         covar.train=data.matrix(ucsf_covs[, c(1:5)])) 

system("rm ~/temp_matrix.bk")
fitdat_test=as_FBM(scale(uq), backingfile = "~/temp_matrix")

ucsf.uq.predictions = predict(res.all.ucsf, fitdat_test, 
                          covar.row=data.matrix(uq_covs[, c(1:5)])) 

print(pROC::auc(uq_status, ucsf.uq.predictions))

########################################################
##### UQ trained UCSF tested

system("rm ~/temp_matrix.bk")
fitdat=as_FBM(scale(uq), backingfile = "~/temp_matrix")

res.all.uq=big_spLogReg(X=fitdat, y01.train=uq_status, alpha=1e-4, 
                    covar.train=data.matrix(uq_covs[, c(1:4)])) ## don't use race since UQ is all white 

system("rm ~/temp_matrix.bk")
fitdat_test=as_FBM(scale(ucsf), backingfile = "~/temp_matrix")

uq.ucsf.predictions = predict(res.all.uq, fitdat_test, 
                           covar.row=data.matrix(ucsf_covs[, c(1:4)])) 

print(pROC::auc(ucsf_status, uq.ucsf.predictions))

########################################################
##### UQ ten fold CV 

system("rm ~/temp_matrix.bk")
fitdat=as_FBM(uq, backingfile = "~/temp_matrix")

nfolds=10
testinds=caret::createFolds(y=uq_status, k=nfolds, list = T)
cvpredicted.uq.only=rep(NA, nfolds)
uq_to_uq <- data.frame(matrix(NA, nrow=nrow(uq), ncol=10))


for(fold in 1:nfolds){
    message(fold)
    traininds=setdiff(1:nrow(uq), testinds[[fold]])
    
    res=big_spLogReg(X=fitdat, y01.train=uq_status[traininds], 
                     alphas=c(1e-4), covar.train=data.matrix(uq_covs[traininds, c(1:4)]),
                    ind.train=traininds, warn=FALSE)
    cvpredicted.uq.only[testinds[[fold]]]=predict(res, fitdat, 
                                                  ind.row=testinds[[fold]], 
                                                 covar.row=data.matrix(uq_covs[testinds[[fold]],
                                                                               c(1:4)]))
}

print(pROC::auc(uq_status, cvpredicted.uq.only))


########################################################
##### UCSF ten fold CV 

system("rm ~/temp_matrix.bk")
fitdat=as_FBM(ucsf, backingfile = "~/temp_matrix")

nfolds=10
testinds=caret::createFolds(y=ucsf_status, k=nfolds, list = T)
cvpredicted.ucsf.only=rep(NA, nfolds)

for(fold in 1:nfolds){
    message(fold)
    traininds=setdiff(1:nrow(ucsf), testinds[[fold]])
    
    res=big_spLogReg(X=fitdat, y01.train=ucsf_status[traininds], 
                     alphas=c(1e-4), covar.train=data.matrix(ucsf_covs[traininds, c(1:5)]),
                    ind.train=traininds, warn=FALSE)
    
    cvpredicted.ucsf.only[testinds[[fold]]]=predict(res, fitdat, ind.row=testinds[[fold]], 
                                        covar.row=data.matrix(ucsf_covs[testinds[[fold]], c(1:5)]))
}
print(pROC::auc(ucsf_status, cvpredicted.ucsf.only))

