# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hao.test<-function(a,miu,cov,alpha){##在本函数中需要的参数是样本数据矩阵，待检验均值(需要是行向量)，已知的协方差阵和设定的显著性水平
  p<-ncol(a)
  mean<-apply(a,2,mean)##a是待检验总体
  mean<-as.matrix(t(mean))
  cov<-as.matrix(cov)
  miu<-as.matrix(miu)
  X2<-4*(mean-miu)%*%cov^(-1)%*%(t(mean-miu))##这里的cov(a)是总体的协方差矩阵
  X2<-as.numeric(X2)
  if(X2>qchisq(1-alpha,p))
    print("均值miu不等于miu0")
  else
    print("均值miu等于miu0")
  print(paste("X2统计量为",X2))
  print(paste("在给定显著性水平下卡方临界值为",qchisq(1-alpha,p)))
}

tang.test<-function(a,miu){
  mean<-apply(a,2,mean)##a是待检验总体
  mean<-as.matrix(t(mean))
  n<-nrow(a)
  L=(n-1)*cov(a)##在这里，n为样本个数
  miu<-as.matrix(miu)
  T2<-n*(n-1)*(mean-miu)%*%(L^(-1))%*%(t(mean-miu))
  print(T2)
}

lv.test<-function(x1,x2,alpha){
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  meanx1<-apply(x1,2,mean)
  meanx1<-t(meanx1)
  meanx1<-as.matrix(meanx1)
  meanx2<-apply(x2,2,mean)
  meanx2<-t(meanx2)
  meanx2<-as.matrix(meanx2)
  cov<-(cov(x1)*(nrow(x1)-1)+cov(x2)*(nrow(x2)-1))/(nrow(x1)+nrow(x2)-2)
  T2<-((nrow(x1)*nrow(x2))/(nrow(x1)+nrow(x2)))*(meanx1-meanx2)%*%cov^(-1)%*%(t(meanx1-meanx2))
  Fx<-((nrow(x1)+nrow(x2)-ncol(x1)-1)/((nrow(x1)+nrow(x2)-2)*ncol(x1)))*T2##两总体均值比较检验量
  Fx<-as.numeric(Fx)
  F<-qf(alpha,ncol(x1), nrow(x1)+nrow(x2)-ncol(x1)-1,lower.tail = FALSE)##alpha是设定的显著性水平，df1=p（即变量个数），df2=n1+n2-p-1)
  ##比较检验统计量Fx和F临界值，若Fx>F临界值，则拒绝原假设，认为两总体均值不相等，反之则没有充分理由拒绝原假设，认为两总体均值相等
  if(Fx>F)
    print("两总体均值不相等")
  else
    print("两总体均值相等")
  print(paste("F统计量为",Fx))
  print(paste("在给定显著性水平下F的临界值为",F))
}

yang.test<-function(x1,x2,alpha){
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  meanx1<-apply(x1,2,mean)
  meanx1<-t(meanx1)
  meanx1<-as.matrix(meanx1)
  meanx2<-apply(x2,2,mean)
  meanx2<-t(meanx2)
  meanx2<-as.matrix(meanx2)
  S<-((nrow(x1)-1)*cov(x1))/(nrow(x1)*(nrow(x1)-1))+((nrow(x2)-1)*cov(x2))/(nrow(x2)*(nrow(x2)-1))
  T2<-(meanx1-meanx2)%*%S^(-1)%*%(t(meanx1-meanx2))##可能出现T^2<0的情况，现在未知如何应对这种情况
  finverse=((nrow(x1)^3-nrow(x1)^2)^(-1))*(((meanx1-meanx2)%*%S^(-1)%*%cov(x1)%*%S^(-1)%*%(t(meanx1-meanx2)))^2)*T2^(-1/8)+((nrow(x2)^3-nrow(x2)^2)^(-1))*(((meanx1-meanx2)%*%S^(-1)%*%cov(x2)%*%S^(-1)%*%(t(meanx1-meanx2)))^2)*T2^(-1/8)
  f<-1/finverse
  Fx<-(f-ncol(x1)+1)/(f*ncol(x1))*T2##检验统计量
  F<-qf(alpha,ncol(x1), f-ncol(x1)+1,lower.tail = FALSE)
  if(Fx>F)
    print("两总体均值不相等")
  else
    print("两总体均值相等")
  print(paste("F统计量为",Fx))
  print(paste("在给定显著性水平下F的临界值为",F))
}

xia.test<-function(a,sigma0,alpha){##sigma0是已知的协方差阵
  sigma0<-as.matrix(sigma0)
  n<-nrow(a)##n是样本个数
  p<-ncol(a)##p是指标个数
  A<-cov(a)%*%sigma0^(-1)
  M<-(n-1)*(ln(det(sigma0))-p-ln(det(cov(a)))+sum(diag(A)))
  D1<-(2*p+1-(2/(p+1)))/(6*(n-1))
  D2<-(p-1)*(p+2)/(6*(n-1)^2)
  f1<-p*(p+1)/2
  f2<-(f1+2)/(D2-D1^2)
  b<-f1/(1-D1-(f1/f2))
  Fx<-M/b
  F<-qf(alpha,f1, f2,lower.tail = FALSE)##df1=f1,df2=f2,显著性水平为alpha为f临界值，F>F临界值，拒绝原假设，认为协方差不相等；反之则没有充分理由拒绝原假设，认为二者协方差相等
  if(Fx>F)
    print("sigma不等于sigma0")
  else
    print("sigma等于sigma0")
  print(paste("F统计量为",Fx))
  print(paste("在给定显著性水平下F的临界值为",F))
}

pcascore<-function(a){
  library(psych)
  a<-scale(a)##对原数据作标准化
  n<-nrow(a)##样本数量
  j<-ncol(a)##属性变量个数
  fa.parallel(a,fa="pc",n.iter=100)##生成碎石图
  pca <- principal(a, rotate = "none")##作主成分分析
  lamda<-data.frame(pca$values)##特征根的值
  p<-sum(rowSums(as.matrix(lamda) > 1))##计算特征根大于1的个数，并将其作为提取主成分的个数
  pca <- principal(a, nfactors = p,rotate = "none")##作主成分分析
  pca##输出结果中，SS loadings中数值表示提取出的两个主成分所对应的特征根大小（在这里分别是6.15和1.47，也可以从上面一行代码的pca$values中找到每一个特征根的值）
  ##Proportion Var表示每一个主成分的方差贡献率大小
  ##Cumulative Var表示主成分的累计方差贡献率
  ##Proportion Explained表示每个主成分在提取的所有主成分中所贡献的方差贡献率大小

  ##在pca输出的数据框中展现的是因子载荷矩阵,即pca$loadings
  pcastructure<-data.frame(pca$Structure[1:j,1:p])
  pcacoefficient<-data.frame()
  for (i in 1:j) {
    for (j in 1:p) {
      pcacoefficient[i,j]<-pcastructure[i,j]/sqrt(lamda[j,1])
    }
  }
  pcacoefficient##主成分系数矩阵(主成分得分系数矩阵)
  score<-data.frame()
  for (i in 1:n) {
    for (j in 1:p) {
      score[i,j]<-sum(a[i,]*pcacoefficient[,j])
    }
  }##计算主成分得分
  print(paste("提取了",p,"个主成分"))
  print("因子载荷矩阵为")
  print(pca$loadings)
  print(score)
}

fcascore<-function(a){
  library(psych)
  correlations<-cor(a)
  n<-nrow(a)
  j<-ncol(a)##属性变量的个数
  scalea<-scale(a)
  fa.parallel(correlations, n.obs = n, fa = "both", n.iter = 100, main = "平行分析碎石图")##通过碎石图判别需要提取几个因子
  ev <- eigen(correlations)##也可以通过特征根的大小来进行判断
  ev$values##查看有几个特征根大小大于1
  p<-sum(rowSums(as.matrix(ev$values) > 1))##计算特征根大于1的个数，并将其作为提取主成分的个数
  ev$vectors##这是特征根对应的特征向量
  fa <- fa(correlations, nfactors = p, rotate = "none",fm="pa")##采取主轴因子法（fm="pa"），提取3个因子(nfactors=p)且不旋转(rotate="none")##在R语言的因子分析中未找到采用主成分法的因子分析，此法待后续版本作补充
  fa
  ##fa结果中的数据框表示的是成分载荷，表示的是观测变量与因子的相关系数
  ##h2表示公因子方差，即因子对每个变量的方法解释度
  lamda<-data.frame(fa$values)##特征根的值
  fastructure<-data.frame(fa$Structure[1:j,1:p])
  facoefficient<-data.frame()
  for (i in 1:j) {
    for (j in 1:p) {
      facoefficient[i,j]<-fastructure[i,j]/lamda[j,1]
    }
  }
  facoefficient##因子得分系数矩阵
  score<-data.frame()
  correlations<-data.frame(correlations)
  for (i in 1:n) {
    for (j in 1:p) {
      score[i,j]<-sum(scalea[i,]*facoefficient[,j])
    }
  }##计算因子得分
  print(score)##各样本所对应三个公因子得分，建议大家采用SPSS作因子分析，这份程序可能计算结果有误！！！
}
