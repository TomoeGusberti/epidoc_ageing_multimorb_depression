4 Multi-group assessment with Cross-Lagged Structural Equation Modelling
analysis
================
Tomoe Gusberti, Dr. Eng
2025-03-20

load data

``` r
load(file='./Data/Data_Imp20_forSEM.RData')#data.mi
GVAr='Sex'
folderName=paste0('MGroup_',GVAr)
```

# Defining models

## model with restriction on ED

``` r
Model_ED<-'
# evolução
#Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
Dep_EpiDoc4~Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm14*MM_EpiDoc2+mm04*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+0*MM_EpiDoc1
Dep_EpiDoc2~dm21*MM_EpiDoc1
MM_EpiDoc4~0*Dep_EpiDoc2+Dep_EpiDoc1 
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~0*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4


# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~ma4*Age
Dep_EpiDoc2~da2*Age
MM_EpiDoc2~ma2*Age
Dep_EpiDoc1~da1*Age
MM_EpiDoc1~ma1*Age

# Escolaridade
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~de2*Education
Dep_EpiDoc1~Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~0*Education
MM_EpiDoc1~Education

'
Model_ED2<-'
# evolução
Dep_EpiDoc4~c(dd14a,dd14b)*Dep_EpiDoc1+c(dd24a,dd24b)*Dep_EpiDoc2
Dep_EpiDoc2~c(dd12a,dd12b)*Dep_EpiDoc1

MM_EpiDoc4~c(mm24a,mm24b)*MM_EpiDoc2+c(mm14a,mm14b)*MM_EpiDoc1
MM_EpiDoc2~c(mm12a,mm12b)*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~c(md24a,md24b)*MM_EpiDoc2+c(md14a,md14b)*MM_EpiDoc1
Dep_EpiDoc2~c(md12a,md12b)*MM_EpiDoc1

MM_EpiDoc4~c(dm24a,dm24b)*Dep_EpiDoc2+c(dm14a,dm14b)*Dep_EpiDoc1
MM_EpiDoc2~c(dm12a,dm12b)*Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age



# educação
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~c(ed2a,ed2b)*Education
Dep_EpiDoc1~c(ed1a,ed1b)*Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~0*Education
MM_EpiDoc1~0*Education


# effect on multimorbidity


# effects on depression =============
##indirect effect age on dep4
ada:=ad1*dd14a+ad2*dd24a+ad1*dd12a*dd24a
adb:=ad1*dd14b+ad2*dd24b+ad1*dd12b*dd24b
amda:=am1a*md14a+am2a*md24a+am1a*mm12a*md24a+ad1*dm12a*md24a
amdb:=am1b*md14b+am2b*md24b+am1b*mm12b*md24b+ad1*dm12b*md24b


adta:=ada+amda
adtb:=adb+amdb

##indirect effect MM1 on dep4
mda:=md12a*dd24a+mm12a*md24a+md14a
mdb:=md12b*dd24b+mm12b*md24b+md14b

## cumulative effect d1 on dep4
dda:=dm12a*md24a+dd12a*dd24a+dd14a
ddb:=dm12b*md24b+dd12b*dd24b+dd14b

## indirect effect of educaiton on dep4
eda:=ed1a*dd12a*dd24a+ed2a*dd24a
edb:=ed1b*dd12b*dd24b+ed2b*dd24b
emda:=ed1a*dm12a*md24a
emdb:=ed1b*dm12b*md24b

edta:=eda+amda
edtb:=edb+emdb

'
```

``` r
resED<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_ED2,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

## reported on fully adjusted

``` r
partable<-parameterEstimates.mi(resED,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.17 (0.031)      0.351 (0)     [0.015:0.325]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2  0.247 (0.002)  0.227 (0.003)     [0.089:0.405]
    ## 3   Dep_EpiDoc2~Dep_EpiDoc1      0.404 (0)      0.373 (0)     [0.266:0.542]
    ## 4     MM_EpiDoc4~MM_EpiDoc2  0.369 (0.016)  0.285 (0.012)     [0.068:0.671]
    ## 5     MM_EpiDoc4~MM_EpiDoc1       0.56 (0)      0.347 (0)     [0.385:0.736]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.102 (0.029)      0.365 (0)      [0.01:0.194]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.494 (0.083)  0.675 (0.028)    [-0.065:1.053]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1  0.077 (0.697) -0.309 (0.152)    [-0.312:0.466]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1   0.23 (0.299)  0.086 (0.727)    [-0.205:0.664]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2  0.035 (0.244)  0.001 (0.955)    [-0.024:0.094]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1  0.004 (0.879)  0.028 (0.118)    [-0.052:0.061]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.018 (0.228)   0.009 (0.46)    [-0.011:0.048]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1  0.549 (0.011)  0.922 (0.003)     [0.126:0.973]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4         0 (NA)         0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age         0 (NA)         0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age   0.01 (0.463)   0.01 (0.463)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age  0.046 (0.002)  0.046 (0.002)     [0.017:0.074]
    ## 19           MM_EpiDoc4~Age  0.018 (0.011)      0.031 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age      0 (0.914)  -0.005 (0.29)    [-0.003:0.004]
    ## 21           MM_EpiDoc1~Age      0.032 (0)      0.049 (0)     [0.023:0.041]
    ## 22    Dep_EpiDoc4~Education         0 (NA)         0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education   0.49 (0.001)   0.435 (0.04)     [0.202:0.777]
    ## 24    Dep_EpiDoc1~Education  0.448 (0.002)      0.776 (0)     [0.161:0.735]
    ## 25     MM_EpiDoc4~Education         0 (NA)         0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education         0 (NA)         0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education         0 (NA)         0 (NA)             [0:0]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4       6.91 (0)       9.37 (0)     [4.942:8.878]
    ## 29 Dep_EpiDoc2~~Dep_EpiDoc2      6.201 (0)     10.759 (0)     [4.487:7.914]
    ## 30   MM_EpiDoc4~~MM_EpiDoc4      0.971 (0)      0.999 (0)     [0.764:1.178]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2      0.189 (0)      0.401 (0)      [0.087:0.29]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1      8.431 (0)     12.642 (0)    [6.156:10.706]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1      0.701 (0)      0.782 (0)     [0.531:0.871]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1      1.461 (0)      1.293 (0)     [1.042:1.879]
    ## 38            Dep_EpiDoc2~1 -0.261 (0.646)  0.513 (0.401)    [-1.383:0.862]
    ## 39             MM_EpiDoc4~1 -0.362 (0.162) -0.559 (0.015)    [-0.869:0.145]
    ## 40             MM_EpiDoc2~1  -0.02 (0.756)  0.165 (0.278)    [-0.144:0.105]
    ## 41            Dep_EpiDoc1~1  -0.657 (0.29) -0.157 (0.819)    [-1.872:0.559]
    ## 42             MM_EpiDoc1~1     -0.785 (0)     -1.335 (0)   [-1.129:-0.441]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1      [0.205:0.496]
    ## 2       [0.08:0.374]
    ## 3      [0.208:0.538]
    ## 4      [0.064:0.506]
    ## 5      [0.176:0.519]
    ## 6      [0.177:0.553]
    ## 7      [0.073:1.277]
    ## 8     [-0.731:0.114]
    ## 9     [-0.399:0.572]
    ## 10    [-0.032:0.034]
    ## 11    [-0.007:0.062]
    ## 12    [-0.015:0.034]
    ## 13     [0.315:1.528]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.017:0.074]
    ## 19     [0.017:0.044]
    ## 20    [-0.013:0.004]
    ## 21     [0.038:0.059]
    ## 22             [0:0]
    ## 23       [0.02:0.85]
    ## 24     [0.389:1.162]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27             [0:0]
    ## 28       [7.3:11.44]
    ## 29    [8.104:13.414]
    ## 30     [0.833:1.166]
    ## 31      [0.28:0.521]
    ## 32    [9.241:16.043]
    ## 33      [0.543:1.02]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.809:1.776]
    ## 38    [-0.692:1.719]
    ## 39    [-1.01:-0.107]
    ## 40    [-0.133:0.464]
    ## 41      [-1.5:1.186]
    ## 42   [-1.728:-0.941]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

``` r
source('./code/parTableToCSV.R')
parTableToCSV(resED,path=paste0('./Results/finalModel'))
```

    ## Joining with `by = join_by(group)`

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `Group2 = as.numeric(substr(label, 2, 2))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

    ##    Group2 label     estimate_p             CI Group
    ## 1      NA     a  0.015 (0.012)  [0.003:0.026]  <NA>
    ## 2      NA     b  0.022 (0.003)  [0.008:0.037]  <NA>
    ## 3      NA    da  0.005 (0.468) [-0.008:0.017]  <NA>
    ## 4      NA    db -0.006 (0.589) [-0.027:0.016]  <NA>
    ## 5      NA    ta  0.019 (0.006)  [0.006:0.033]  <NA>
    ## 6      NA    tb   0.016 (0.19) [-0.008:0.041]  <NA>
    ## 7      NA     a  0.183 (0.386)  [-0.23:0.596]  <NA>
    ## 8      NA     b -0.044 (0.844) [-0.487:0.399]  <NA>
    ## 9      NA     a      0.279 (0)  [0.125:0.433]  <NA>
    ## 10     NA     b      0.442 (0)  [0.311:0.572]  <NA>
    ## 11     NA     a  0.166 (0.012)  [0.036:0.296]  <NA>
    ## 12     NA     b  0.164 (0.039)  [0.008:0.319]  <NA>
    ## 13     NA    da  0.004 (0.397) [-0.005:0.013]  <NA>
    ## 14     NA    db  0.005 (0.524)   [-0.01:0.02]  <NA>
    ## 15     NA    ta  0.171 (0.009)  [0.042:0.299]  <NA>
    ## 16     NA    tb  0.168 (0.037)   [0.01:0.327]  <NA>

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : there are records with missing times, which will be dropped.

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : multiple rows match for Group=NA: first taken

    ## $partable
    ##             lhs op                                                    rhs group
    ## 1   Dep_EpiDoc4  ~                                            Dep_EpiDoc1     1
    ## 2   Dep_EpiDoc4  ~                                            Dep_EpiDoc2     1
    ## 3   Dep_EpiDoc2  ~                                            Dep_EpiDoc1     1
    ## 4    MM_EpiDoc4  ~                                             MM_EpiDoc2     1
    ## 5    MM_EpiDoc4  ~                                             MM_EpiDoc1     1
    ## 6    MM_EpiDoc2  ~                                             MM_EpiDoc1     1
    ## 7   Dep_EpiDoc4  ~                                             MM_EpiDoc2     1
    ## 8   Dep_EpiDoc4  ~                                             MM_EpiDoc1     1
    ## 9   Dep_EpiDoc2  ~                                             MM_EpiDoc1     1
    ## 10   MM_EpiDoc4  ~                                            Dep_EpiDoc2     1
    ## 11   MM_EpiDoc4  ~                                            Dep_EpiDoc1     1
    ## 12   MM_EpiDoc2  ~                                            Dep_EpiDoc1     1
    ## 13  Dep_EpiDoc1 ~~                                             MM_EpiDoc1     1
    ## 14  Dep_EpiDoc2 ~~                                             MM_EpiDoc2     1
    ## 15  Dep_EpiDoc4 ~~                                             MM_EpiDoc4     1
    ## 16  Dep_EpiDoc4  ~                                                    Age     1
    ## 17  Dep_EpiDoc2  ~                                                    Age     1
    ## 18  Dep_EpiDoc1  ~                                                    Age     1
    ## 19   MM_EpiDoc4  ~                                                    Age     1
    ## 20   MM_EpiDoc2  ~                                                    Age     1
    ## 21   MM_EpiDoc1  ~                                                    Age     1
    ## 22  Dep_EpiDoc4  ~                                              Education     1
    ## 23  Dep_EpiDoc2  ~                                              Education     1
    ## 24  Dep_EpiDoc1  ~                                              Education     1
    ## 25   MM_EpiDoc4  ~                                              Education     1
    ## 26   MM_EpiDoc2  ~                                              Education     1
    ## 27   MM_EpiDoc1  ~                                              Education     1
    ## 28  Dep_EpiDoc4 ~~                                            Dep_EpiDoc4     1
    ## 29  Dep_EpiDoc2 ~~                                            Dep_EpiDoc2     1
    ## 30   MM_EpiDoc4 ~~                                             MM_EpiDoc4     1
    ## 31   MM_EpiDoc2 ~~                                             MM_EpiDoc2     1
    ## 32  Dep_EpiDoc1 ~~                                            Dep_EpiDoc1     1
    ## 33   MM_EpiDoc1 ~~                                             MM_EpiDoc1     1
    ## 34          Age ~~                                                    Age     1
    ## 35          Age ~~                                              Education     1
    ## 36    Education ~~                                              Education     1
    ## 37  Dep_EpiDoc4 ~1                                                            1
    ## 38  Dep_EpiDoc2 ~1                                                            1
    ## 39   MM_EpiDoc4 ~1                                                            1
    ## 40   MM_EpiDoc2 ~1                                                            1
    ## 41  Dep_EpiDoc1 ~1                                                            1
    ## 42   MM_EpiDoc1 ~1                                                            1
    ## 43          Age ~1                                                            1
    ## 44    Education ~1                                                            1
    ## 45  Dep_EpiDoc4  ~                                            Dep_EpiDoc1     2
    ## 46  Dep_EpiDoc4  ~                                            Dep_EpiDoc2     2
    ## 47  Dep_EpiDoc2  ~                                            Dep_EpiDoc1     2
    ## 48   MM_EpiDoc4  ~                                             MM_EpiDoc2     2
    ## 49   MM_EpiDoc4  ~                                             MM_EpiDoc1     2
    ## 50   MM_EpiDoc2  ~                                             MM_EpiDoc1     2
    ## 51  Dep_EpiDoc4  ~                                             MM_EpiDoc2     2
    ## 52  Dep_EpiDoc4  ~                                             MM_EpiDoc1     2
    ## 53  Dep_EpiDoc2  ~                                             MM_EpiDoc1     2
    ## 54   MM_EpiDoc4  ~                                            Dep_EpiDoc2     2
    ## 55   MM_EpiDoc4  ~                                            Dep_EpiDoc1     2
    ## 56   MM_EpiDoc2  ~                                            Dep_EpiDoc1     2
    ## 57  Dep_EpiDoc1 ~~                                             MM_EpiDoc1     2
    ## 58  Dep_EpiDoc2 ~~                                             MM_EpiDoc2     2
    ## 59  Dep_EpiDoc4 ~~                                             MM_EpiDoc4     2
    ## 60  Dep_EpiDoc4  ~                                                    Age     2
    ## 61  Dep_EpiDoc2  ~                                                    Age     2
    ## 62  Dep_EpiDoc1  ~                                                    Age     2
    ## 63   MM_EpiDoc4  ~                                                    Age     2
    ## 64   MM_EpiDoc2  ~                                                    Age     2
    ## 65   MM_EpiDoc1  ~                                                    Age     2
    ## 66  Dep_EpiDoc4  ~                                              Education     2
    ## 67  Dep_EpiDoc2  ~                                              Education     2
    ## 68  Dep_EpiDoc1  ~                                              Education     2
    ## 69   MM_EpiDoc4  ~                                              Education     2
    ## 70   MM_EpiDoc2  ~                                              Education     2
    ## 71   MM_EpiDoc1  ~                                              Education     2
    ## 72  Dep_EpiDoc4 ~~                                            Dep_EpiDoc4     2
    ## 73  Dep_EpiDoc2 ~~                                            Dep_EpiDoc2     2
    ## 74   MM_EpiDoc4 ~~                                             MM_EpiDoc4     2
    ## 75   MM_EpiDoc2 ~~                                             MM_EpiDoc2     2
    ## 76  Dep_EpiDoc1 ~~                                            Dep_EpiDoc1     2
    ## 77   MM_EpiDoc1 ~~                                             MM_EpiDoc1     2
    ## 78          Age ~~                                                    Age     2
    ## 79          Age ~~                                              Education     2
    ## 80    Education ~~                                              Education     2
    ## 81  Dep_EpiDoc4 ~1                                                            2
    ## 82  Dep_EpiDoc2 ~1                                                            2
    ## 83   MM_EpiDoc4 ~1                                                            2
    ## 84   MM_EpiDoc2 ~1                                                            2
    ## 85  Dep_EpiDoc1 ~1                                                            2
    ## 86   MM_EpiDoc1 ~1                                                            2
    ## 87          Age ~1                                                            2
    ## 88    Education ~1                                                            2
    ## 89          ada :=                    ad1*dd14a+ad2*dd24a+ad1*dd12a*dd24a     0
    ## 90          adb :=                    ad1*dd14b+ad2*dd24b+ad1*dd12b*dd24b     0
    ## 91         amda := am1a*md14a+am2a*md24a+am1a*mm12a*md24a+ad1*dm12a*md24a     0
    ## 92         amdb := am1b*md14b+am2b*md24b+am1b*mm12b*md24b+ad1*dm12b*md24b     0
    ## 93         adta :=                                               ada+amda     0
    ## 94         adtb :=                                               adb+amdb     0
    ## 95          mda :=                          md12a*dd24a+mm12a*md24a+md14a     0
    ## 96          mdb :=                          md12b*dd24b+mm12b*md24b+md14b     0
    ## 97          dda :=                          dm12a*md24a+dd12a*dd24a+dd14a     0
    ## 98          ddb :=                          dm12b*md24b+dd12b*dd24b+dd14b     0
    ## 99          eda :=                            ed1a*dd12a*dd24a+ed2a*dd24a     0
    ## 100         edb :=                            ed1b*dd12b*dd24b+ed2b*dd24b     0
    ## 101        emda :=                                       ed1a*dm12a*md24a     0
    ## 102        emdb :=                                       ed1b*dm12b*md24b     0
    ## 103        edta :=                                               eda+amda     0
    ## 104        edtb :=                                               edb+emdb     0
    ##     label block     est    se      t           df pvalue ci.lower ci.upper
    ## 1   dd14a     1   0.170 0.079  2.158 5.639800e+02  0.031    0.015    0.325
    ## 2   dd24a     1   0.247 0.080  3.079 3.345720e+02  0.002    0.089    0.405
    ## 3   dd12a     1   0.404 0.070  5.772 3.829560e+02  0.000    0.266    0.542
    ## 4   mm24a     1   0.369 0.154  2.398 6.993126e+06  0.016    0.068    0.671
    ## 5   mm14a     1   0.560 0.089  6.265 3.170168e+05  0.000    0.385    0.736
    ## 6   mm12a     1   0.102 0.047  2.183 1.124620e+25  0.029    0.010    0.194
    ## 7   md24a     1   0.494 0.285  1.732 2.055989e+03  0.083   -0.065    1.053
    ## 8   md14a     1   0.077 0.198  0.389 1.081051e+03  0.697   -0.312    0.466
    ## 9   md12a     1   0.230 0.221  1.039 3.582280e+02  0.299   -0.205    0.664
    ## 10  dm24a     1   0.035 0.030  1.169 2.104400e+02  0.244   -0.024    0.094
    ## 11  dm14a     1   0.004 0.029  0.152 2.472136e+03  0.879   -0.052    0.061
    ## 12  dm12a     1   0.018 0.015  1.206 2.618237e+27  0.228   -0.011    0.048
    ## 13            1   0.549 0.216  2.543 4.011747e+23  0.011    0.126    0.973
    ## 14            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 15            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 16            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 17    ad2     1   0.010 0.013  0.737 8.894000e+01  0.463   -0.016    0.035
    ## 18    ad1     1   0.046 0.014  3.165 5.335687e+23  0.002    0.017    0.074
    ## 19   am4a     1   0.018 0.007  2.534 1.489678e+07  0.011    0.004    0.033
    ## 20   am2a     1   0.000 0.002  0.108 3.012649e+22  0.914   -0.003    0.004
    ## 21   am1a     1   0.032 0.005  7.009 1.655413e+21  0.000    0.023    0.041
    ## 22            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 23   ed2a     1   0.490 0.146  3.345 3.353830e+02  0.001    0.202    0.777
    ## 24   ed1a     1   0.448 0.146  3.061 7.048221e+22  0.002    0.161    0.735
    ## 25            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 26            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 27            1   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 28            1   6.910 1.003  6.887 1.748822e+03  0.000    4.942    8.878
    ## 29            1   6.201 0.864  7.178 1.032680e+02  0.000    4.487    7.914
    ## 30            1   0.971 0.106  9.179 5.893331e+05  0.000    0.764    1.178
    ## 31            1   0.189 0.052  3.638 5.210387e+28  0.000    0.087    0.290
    ## 32            1   8.431 1.161  7.264 1.705387e+24  0.000    6.156   10.706
    ## 33            1   0.701 0.087  8.072 1.903381e+23  0.000    0.531    0.871
    ## 34            1 204.021 0.000     NA           NA     NA  204.021  204.021
    ## 35            1   6.794 0.000     NA           NA     NA    6.794    6.794
    ## 36            1   1.159 0.000     NA           NA     NA    1.159    1.159
    ## 37            1   1.461 0.213  6.859 4.014490e+02  0.000    1.042    1.879
    ## 38            1  -0.261 0.567 -0.460 1.225420e+02  0.646   -1.383    0.862
    ## 39            1  -0.362 0.259 -1.400 2.410734e+05  0.162   -0.869    0.145
    ## 40            1  -0.020 0.063 -0.311 1.124786e+22  0.756   -0.144    0.105
    ## 41            1  -0.657 0.620 -1.059 7.143187e+24  0.290   -1.872    0.559
    ## 42            1  -0.785 0.176 -4.469 6.013346e+20  0.000   -1.129   -0.441
    ## 43            1  44.188 0.000     NA           NA     NA   44.188   44.188
    ## 44            1   2.823 0.000     NA           NA     NA    2.823    2.823
    ## 45  dd14b     2   0.351 0.074  4.710 1.112120e+04  0.000    0.205    0.496
    ## 46  dd24b     2   0.227 0.075  3.024 1.157530e+03  0.003    0.080    0.374
    ## 47  dd12b     2   0.373 0.084  4.444 8.443920e+02  0.000    0.208    0.538
    ## 48  mm24b     2   0.285 0.113  2.524 3.404261e+07  0.012    0.064    0.506
    ## 49  mm14b     2   0.347 0.088  3.968 2.778892e+08  0.000    0.176    0.519
    ## 50  mm12b     2   0.365 0.096  3.803 4.422834e+23  0.000    0.177    0.553
    ## 51  md24b     2   0.675 0.307  2.199 1.816978e+03  0.028    0.073    1.277
    ## 52  md14b     2  -0.309 0.215 -1.433 1.722086e+03  0.152   -0.731    0.114
    ## 53  md12b     2   0.086 0.247  0.349 3.012690e+02  0.727   -0.399    0.572
    ## 54  dm24b     2   0.001 0.017  0.057 4.923450e+02  0.955   -0.032    0.034
    ## 55  dm14b     2   0.028 0.018  1.562 1.503371e+04  0.118   -0.007    0.062
    ## 56  dm12b     2   0.009 0.013  0.739 1.882637e+24  0.460   -0.015    0.034
    ## 57            2   0.922 0.309  2.979 7.852117e+24  0.003    0.315    1.528
    ## 58            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 59            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 60            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 61    ad2     2   0.010 0.013  0.737 8.894000e+01  0.463   -0.016    0.035
    ## 62    ad1     2   0.046 0.014  3.165 5.335687e+23  0.002    0.017    0.074
    ## 63   am4b     2   0.031 0.007  4.502 9.888082e+05  0.000    0.017    0.044
    ## 64   am2b     2  -0.005 0.004 -1.059 1.586137e+21  0.290   -0.013    0.004
    ## 65   am1b     2   0.049 0.005  9.116 1.060970e+21  0.000    0.038    0.059
    ## 66            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 67   ed2b     2   0.435 0.211  2.061 2.963950e+02  0.040    0.020    0.850
    ## 68   ed1b     2   0.776 0.197  3.935 6.479977e+22  0.000    0.389    1.162
    ## 69            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 70            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 71            2   0.000 0.000     NA           NA     NA    0.000    0.000
    ## 72            2   9.370 1.055  8.880 1.277156e+03  0.000    7.300   11.440
    ## 73            2  10.759 1.350  7.968 3.790090e+02  0.000    8.104   13.414
    ## 74            2   0.999 0.085 11.762 1.131643e+09  0.000    0.833    1.166
    ## 75            2   0.401 0.061  6.524 1.238205e+30  0.000    0.280    0.521
    ## 76            2  12.642 1.735  7.286 1.892209e+24  0.000    9.241   16.043
    ## 77            2   0.782 0.122  6.411 3.233451e+23  0.000    0.543    1.020
    ## 78            2 135.726 0.000     NA           NA     NA  135.726  135.726
    ## 79            2   3.680 0.000     NA           NA     NA    3.680    3.680
    ## 80            2   1.240 0.000     NA           NA     NA    1.240    1.240
    ## 81            2   1.293 0.247  5.244 3.903109e+03  0.000    0.809    1.776
    ## 82            2   0.513 0.610  0.841 1.470740e+02  0.401   -0.692    1.719
    ## 83            2  -0.559 0.230 -2.425 1.438792e+08  0.015   -1.010   -0.107
    ## 84            2   0.165 0.152  1.084 6.420069e+20  0.278   -0.133    0.464
    ## 85            2  -0.157 0.685 -0.229 8.262257e+23  0.819   -1.500    1.186
    ## 86            2  -1.335 0.201 -6.648 5.044878e+20  0.000   -1.728   -0.941
    ## 87            2  42.481 0.000     NA           NA     NA   42.481   42.481
    ## 88            2   2.409 0.000     NA           NA     NA    2.409    2.409
    ## 89    ada     0   0.015 0.006  2.518 4.171710e+02  0.012    0.003    0.026
    ## 90    adb     0   0.022 0.007  3.027 1.402441e+03  0.003    0.008    0.037
    ## 91   amda     0   0.005 0.006  0.725 1.402271e+03  0.468   -0.008    0.017
    ## 92   amdb     0  -0.006 0.011 -0.541 3.486222e+03  0.589   -0.027    0.016
    ## 93   adta     0   0.019 0.007  2.784 4.396330e+02  0.006    0.006    0.033
    ## 94   adtb     0   0.016 0.012  1.311 1.868028e+03  0.190   -0.008    0.041
    ## 95    mda     0   0.183 0.211  0.867 2.557159e+03  0.386   -0.230    0.596
    ## 96    mdb     0  -0.044 0.226 -0.196 3.119806e+03  0.844   -0.487    0.399
    ## 97    dda     0   0.279 0.079  3.551 2.254727e+03  0.000    0.125    0.433
    ## 98    ddb     0   0.442 0.067  6.625 1.072315e+05  0.000    0.311    0.572
    ## 99    eda     0   0.166 0.066  2.523 2.875870e+02  0.012    0.036    0.296
    ## 100   edb     0   0.164 0.079  2.066 1.909250e+03  0.039    0.008    0.319
    ## 101  emda     0   0.004 0.005  0.846 2.540594e+04  0.397   -0.005    0.013
    ## 102  emdb     0   0.005 0.008  0.637 1.456509e+05  0.524   -0.010    0.020
    ## 103  edta     0   0.171 0.065  2.614 2.902320e+02  0.009    0.042    0.299
    ## 104  edtb     0   0.168 0.081  2.086 2.044522e+03  0.037    0.010    0.327
    ##      std.lv std.all std.nox
    ## 1     0.170   0.181   0.181
    ## 2     0.247   0.255   0.255
    ## 3     0.404   0.417   0.417
    ## 4     0.369   0.131   0.131
    ## 5     0.560   0.421   0.421
    ## 6     0.102   0.216   0.216
    ## 7     0.494   0.078   0.078
    ## 8     0.077   0.026   0.026
    ## 9     0.230   0.074   0.074
    ## 10    0.035   0.082   0.082
    ## 11    0.004   0.011   0.011
    ## 12    0.018   0.123   0.123
    ## 13    0.549   0.226   0.226
    ## 14    0.000   0.000   0.000
    ## 15    0.000   0.000   0.000
    ## 16    0.000   0.000   0.000
    ## 17    0.010   0.046   0.003
    ## 18    0.046   0.214   0.015
    ## 19    0.018   0.206   0.014
    ## 20    0.000   0.006   0.000
    ## 21    0.032   0.483   0.034
    ## 22    0.000   0.000   0.000
    ## 23    0.490   0.177   0.165
    ## 24    0.448   0.158   0.146
    ## 25    0.000   0.000   0.000
    ## 26    0.000   0.000   0.000
    ## 27    0.000   0.000   0.000
    ## 28    6.910   0.832   0.832
    ## 29    6.201   0.703   0.703
    ## 30    0.971   0.599   0.599
    ## 31    0.189   0.919   0.919
    ## 32    8.431   0.900   0.900
    ## 33    0.701   0.767   0.767
    ## 34  204.021   1.000 204.021
    ## 35    6.794   0.442   6.794
    ## 36    1.159   1.000   1.159
    ## 37    1.461   0.507   0.507
    ## 38   -0.261  -0.088  -0.088
    ## 39   -0.362  -0.285  -0.285
    ## 40   -0.020  -0.044  -0.044
    ## 41   -0.657  -0.214  -0.214
    ## 42   -0.785  -0.821  -0.821
    ## 43   44.188   3.094  44.188
    ## 44    2.823   2.622   2.823
    ## 45    0.351   0.365   0.365
    ## 46    0.227   0.233   0.233
    ## 47    0.373   0.379   0.379
    ## 48    0.285   0.164   0.164
    ## 49    0.347   0.287   0.287
    ## 50    0.365   0.522   0.522
    ## 51    0.675   0.138   0.138
    ## 52   -0.309  -0.090  -0.090
    ## 53    0.086   0.025   0.025
    ## 54    0.001   0.003   0.003
    ## 55    0.028   0.081   0.081
    ## 56    0.009   0.047   0.047
    ## 57    0.922   0.293   0.293
    ## 58    0.000   0.000   0.000
    ## 59    0.000   0.000   0.000
    ## 60    0.000   0.000   0.000
    ## 61    0.010   0.030   0.003
    ## 62    0.046   0.143   0.012
    ## 63    0.031   0.280   0.024
    ## 64   -0.005  -0.075  -0.006
    ## 65    0.049   0.539   0.046
    ## 66    0.000   0.000   0.000
    ## 67    0.435   0.132   0.118
    ## 68    0.776   0.231   0.208
    ## 69    0.000   0.000   0.000
    ## 70    0.000   0.000   0.000
    ## 71    0.000   0.000   0.000
    ## 72    9.370   0.728   0.728
    ## 73   10.759   0.795   0.795
    ## 74    0.999   0.620   0.620
    ## 75    0.401   0.746   0.746
    ## 76   12.642   0.907   0.907
    ## 77    0.782   0.709   0.709
    ## 78  135.726   1.000 135.726
    ## 79    3.680   0.284   3.680
    ## 80    1.240   1.000   1.240
    ## 81    1.293   0.360   0.360
    ## 82    0.513   0.140   0.140
    ## 83   -0.559  -0.440  -0.440
    ## 84    0.165   0.225   0.225
    ## 85   -0.157  -0.042  -0.042
    ## 86   -1.335  -1.272  -1.272
    ## 87   42.481   3.646  42.481
    ## 88    2.409   2.163   2.409
    ## 89    0.015   0.073   0.005
    ## 90    0.022   0.108   0.008
    ## 91    0.005   0.023   0.002
    ## 92   -0.006  -0.019  -0.002
    ## 93    0.019   0.096   0.007
    ## 94    0.016   0.089   0.006
    ## 95    0.184   0.061   0.061
    ## 96   -0.043  -0.013  -0.013
    ## 97    0.279   0.296   0.296
    ## 98    0.442   0.459   0.459
    ## 99    0.166   0.062   0.058
    ## 100   0.165   0.051   0.046
    ## 101   0.004   0.002   0.001
    ## 102   0.005   0.002   0.001
    ## 103   0.171   0.085   0.059
    ## 104   0.169   0.053   0.047
    ##                                                           label2     estimate_p
    ## 1                                        Dep_EpiDoc4~Dep_EpiDoc1   0.17 (0.031)
    ## 2                                        Dep_EpiDoc4~Dep_EpiDoc2  0.247 (0.002)
    ## 3                                        Dep_EpiDoc2~Dep_EpiDoc1      0.404 (0)
    ## 4                                          MM_EpiDoc4~MM_EpiDoc2  0.369 (0.016)
    ## 5                                          MM_EpiDoc4~MM_EpiDoc1       0.56 (0)
    ## 6                                          MM_EpiDoc2~MM_EpiDoc1  0.102 (0.029)
    ## 7                                         Dep_EpiDoc4~MM_EpiDoc2  0.494 (0.083)
    ## 8                                         Dep_EpiDoc4~MM_EpiDoc1  0.077 (0.697)
    ## 9                                         Dep_EpiDoc2~MM_EpiDoc1   0.23 (0.299)
    ## 10                                        MM_EpiDoc4~Dep_EpiDoc2  0.035 (0.244)
    ## 11                                        MM_EpiDoc4~Dep_EpiDoc1  0.004 (0.879)
    ## 12                                        MM_EpiDoc2~Dep_EpiDoc1  0.018 (0.228)
    ## 13                                       Dep_EpiDoc1~~MM_EpiDoc1  0.549 (0.011)
    ## 14                                       Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)
    ## 15                                       Dep_EpiDoc4~~MM_EpiDoc4         0 (NA)
    ## 16                                               Dep_EpiDoc4~Age         0 (NA)
    ## 17                                               Dep_EpiDoc2~Age   0.01 (0.463)
    ## 18                                               Dep_EpiDoc1~Age  0.046 (0.002)
    ## 19                                                MM_EpiDoc4~Age  0.018 (0.011)
    ## 20                                                MM_EpiDoc2~Age      0 (0.914)
    ## 21                                                MM_EpiDoc1~Age      0.032 (0)
    ## 22                                         Dep_EpiDoc4~Education         0 (NA)
    ## 23                                         Dep_EpiDoc2~Education   0.49 (0.001)
    ## 24                                         Dep_EpiDoc1~Education  0.448 (0.002)
    ## 25                                          MM_EpiDoc4~Education         0 (NA)
    ## 26                                          MM_EpiDoc2~Education         0 (NA)
    ## 27                                          MM_EpiDoc1~Education         0 (NA)
    ## 28                                      Dep_EpiDoc4~~Dep_EpiDoc4       6.91 (0)
    ## 29                                      Dep_EpiDoc2~~Dep_EpiDoc2      6.201 (0)
    ## 30                                        MM_EpiDoc4~~MM_EpiDoc4      0.971 (0)
    ## 31                                        MM_EpiDoc2~~MM_EpiDoc2      0.189 (0)
    ## 32                                      Dep_EpiDoc1~~Dep_EpiDoc1      8.431 (0)
    ## 33                                        MM_EpiDoc1~~MM_EpiDoc1      0.701 (0)
    ## 34                                                      Age~~Age   204.021 (NA)
    ## 35                                                Age~~Education     6.794 (NA)
    ## 36                                          Education~~Education     1.159 (NA)
    ## 37                                                 Dep_EpiDoc4~1      1.461 (0)
    ## 38                                                 Dep_EpiDoc2~1 -0.261 (0.646)
    ## 39                                                  MM_EpiDoc4~1 -0.362 (0.162)
    ## 40                                                  MM_EpiDoc2~1  -0.02 (0.756)
    ## 41                                                 Dep_EpiDoc1~1  -0.657 (0.29)
    ## 42                                                  MM_EpiDoc1~1     -0.785 (0)
    ## 43                                                         Age~1    44.188 (NA)
    ## 44                                                   Education~1     2.823 (NA)
    ## 45                                       Dep_EpiDoc4~Dep_EpiDoc1      0.351 (0)
    ## 46                                       Dep_EpiDoc4~Dep_EpiDoc2  0.227 (0.003)
    ## 47                                       Dep_EpiDoc2~Dep_EpiDoc1      0.373 (0)
    ## 48                                         MM_EpiDoc4~MM_EpiDoc2  0.285 (0.012)
    ## 49                                         MM_EpiDoc4~MM_EpiDoc1      0.347 (0)
    ## 50                                         MM_EpiDoc2~MM_EpiDoc1      0.365 (0)
    ## 51                                        Dep_EpiDoc4~MM_EpiDoc2  0.675 (0.028)
    ## 52                                        Dep_EpiDoc4~MM_EpiDoc1 -0.309 (0.152)
    ## 53                                        Dep_EpiDoc2~MM_EpiDoc1  0.086 (0.727)
    ## 54                                        MM_EpiDoc4~Dep_EpiDoc2  0.001 (0.955)
    ## 55                                        MM_EpiDoc4~Dep_EpiDoc1  0.028 (0.118)
    ## 56                                        MM_EpiDoc2~Dep_EpiDoc1   0.009 (0.46)
    ## 57                                       Dep_EpiDoc1~~MM_EpiDoc1  0.922 (0.003)
    ## 58                                       Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)
    ## 59                                       Dep_EpiDoc4~~MM_EpiDoc4         0 (NA)
    ## 60                                               Dep_EpiDoc4~Age         0 (NA)
    ## 61                                               Dep_EpiDoc2~Age   0.01 (0.463)
    ## 62                                               Dep_EpiDoc1~Age  0.046 (0.002)
    ## 63                                                MM_EpiDoc4~Age      0.031 (0)
    ## 64                                                MM_EpiDoc2~Age  -0.005 (0.29)
    ## 65                                                MM_EpiDoc1~Age      0.049 (0)
    ## 66                                         Dep_EpiDoc4~Education         0 (NA)
    ## 67                                         Dep_EpiDoc2~Education   0.435 (0.04)
    ## 68                                         Dep_EpiDoc1~Education      0.776 (0)
    ## 69                                          MM_EpiDoc4~Education         0 (NA)
    ## 70                                          MM_EpiDoc2~Education         0 (NA)
    ## 71                                          MM_EpiDoc1~Education         0 (NA)
    ## 72                                      Dep_EpiDoc4~~Dep_EpiDoc4       9.37 (0)
    ## 73                                      Dep_EpiDoc2~~Dep_EpiDoc2     10.759 (0)
    ## 74                                        MM_EpiDoc4~~MM_EpiDoc4      0.999 (0)
    ## 75                                        MM_EpiDoc2~~MM_EpiDoc2      0.401 (0)
    ## 76                                      Dep_EpiDoc1~~Dep_EpiDoc1     12.642 (0)
    ## 77                                        MM_EpiDoc1~~MM_EpiDoc1      0.782 (0)
    ## 78                                                      Age~~Age   135.726 (NA)
    ## 79                                                Age~~Education      3.68 (NA)
    ## 80                                          Education~~Education      1.24 (NA)
    ## 81                                                 Dep_EpiDoc4~1      1.293 (0)
    ## 82                                                 Dep_EpiDoc2~1  0.513 (0.401)
    ## 83                                                  MM_EpiDoc4~1 -0.559 (0.015)
    ## 84                                                  MM_EpiDoc2~1  0.165 (0.278)
    ## 85                                                 Dep_EpiDoc1~1 -0.157 (0.819)
    ## 86                                                  MM_EpiDoc1~1     -1.335 (0)
    ## 87                                                         Age~1    42.481 (NA)
    ## 88                                                   Education~1     2.409 (NA)
    ## 89                      ada:=ad1*dd14a+ad2*dd24a+ad1*dd12a*dd24a  0.015 (0.012)
    ## 90                      adb:=ad1*dd14b+ad2*dd24b+ad1*dd12b*dd24b  0.022 (0.003)
    ## 91  amda:=am1a*md14a+am2a*md24a+am1a*mm12a*md24a+ad1*dm12a*md24a  0.005 (0.468)
    ## 92  amdb:=am1b*md14b+am2b*md24b+am1b*mm12b*md24b+ad1*dm12b*md24b -0.006 (0.589)
    ## 93                                                adta:=ada+amda  0.019 (0.006)
    ## 94                                                adtb:=adb+amdb   0.016 (0.19)
    ## 95                            mda:=md12a*dd24a+mm12a*md24a+md14a  0.183 (0.386)
    ## 96                            mdb:=md12b*dd24b+mm12b*md24b+md14b -0.044 (0.844)
    ## 97                            dda:=dm12a*md24a+dd12a*dd24a+dd14a      0.279 (0)
    ## 98                            ddb:=dm12b*md24b+dd12b*dd24b+dd14b      0.442 (0)
    ## 99                              eda:=ed1a*dd12a*dd24a+ed2a*dd24a  0.166 (0.012)
    ## 100                             edb:=ed1b*dd12b*dd24b+ed2b*dd24b  0.164 (0.039)
    ## 101                                       emda:=ed1a*dm12a*md24a  0.004 (0.397)
    ## 102                                       emdb:=ed1b*dm12b*md24b  0.005 (0.524)
    ## 103                                               edta:=eda+amda  0.171 (0.009)
    ## 104                                               edtb:=edb+emdb  0.168 (0.037)
    ##                    CI  Group
    ## 1       [0.015:0.325]   Male
    ## 2       [0.089:0.405]   Male
    ## 3       [0.266:0.542]   Male
    ## 4       [0.068:0.671]   Male
    ## 5       [0.385:0.736]   Male
    ## 6        [0.01:0.194]   Male
    ## 7      [-0.065:1.053]   Male
    ## 8      [-0.312:0.466]   Male
    ## 9      [-0.205:0.664]   Male
    ## 10     [-0.024:0.094]   Male
    ## 11     [-0.052:0.061]   Male
    ## 12     [-0.011:0.048]   Male
    ## 13      [0.126:0.973]   Male
    ## 14              [0:0]   Male
    ## 15              [0:0]   Male
    ## 16              [0:0]   Male
    ## 17     [-0.016:0.035]   Male
    ## 18      [0.017:0.074]   Male
    ## 19      [0.004:0.033]   Male
    ## 20     [-0.003:0.004]   Male
    ## 21      [0.023:0.041]   Male
    ## 22              [0:0]   Male
    ## 23      [0.202:0.777]   Male
    ## 24      [0.161:0.735]   Male
    ## 25              [0:0]   Male
    ## 26              [0:0]   Male
    ## 27              [0:0]   Male
    ## 28      [4.942:8.878]   Male
    ## 29      [4.487:7.914]   Male
    ## 30      [0.764:1.178]   Male
    ## 31       [0.087:0.29]   Male
    ## 32     [6.156:10.706]   Male
    ## 33      [0.531:0.871]   Male
    ## 34  [204.021:204.021]   Male
    ## 35      [6.794:6.794]   Male
    ## 36      [1.159:1.159]   Male
    ## 37      [1.042:1.879]   Male
    ## 38     [-1.383:0.862]   Male
    ## 39     [-0.869:0.145]   Male
    ## 40     [-0.144:0.105]   Male
    ## 41     [-1.872:0.559]   Male
    ## 42    [-1.129:-0.441]   Male
    ## 43    [44.188:44.188]   Male
    ## 44      [2.823:2.823]   Male
    ## 45      [0.205:0.496] Female
    ## 46       [0.08:0.374] Female
    ## 47      [0.208:0.538] Female
    ## 48      [0.064:0.506] Female
    ## 49      [0.176:0.519] Female
    ## 50      [0.177:0.553] Female
    ## 51      [0.073:1.277] Female
    ## 52     [-0.731:0.114] Female
    ## 53     [-0.399:0.572] Female
    ## 54     [-0.032:0.034] Female
    ## 55     [-0.007:0.062] Female
    ## 56     [-0.015:0.034] Female
    ## 57      [0.315:1.528] Female
    ## 58              [0:0] Female
    ## 59              [0:0] Female
    ## 60              [0:0] Female
    ## 61     [-0.016:0.035] Female
    ## 62      [0.017:0.074] Female
    ## 63      [0.017:0.044] Female
    ## 64     [-0.013:0.004] Female
    ## 65      [0.038:0.059] Female
    ## 66              [0:0] Female
    ## 67        [0.02:0.85] Female
    ## 68      [0.389:1.162] Female
    ## 69              [0:0] Female
    ## 70              [0:0] Female
    ## 71              [0:0] Female
    ## 72        [7.3:11.44] Female
    ## 73     [8.104:13.414] Female
    ## 74      [0.833:1.166] Female
    ## 75       [0.28:0.521] Female
    ## 76     [9.241:16.043] Female
    ## 77       [0.543:1.02] Female
    ## 78  [135.726:135.726] Female
    ## 79        [3.68:3.68] Female
    ## 80        [1.24:1.24] Female
    ## 81      [0.809:1.776] Female
    ## 82     [-0.692:1.719] Female
    ## 83     [-1.01:-0.107] Female
    ## 84     [-0.133:0.464] Female
    ## 85       [-1.5:1.186] Female
    ## 86    [-1.728:-0.941] Female
    ## 87    [42.481:42.481] Female
    ## 88      [2.409:2.409] Female
    ## 89      [0.003:0.026]   <NA>
    ## 90      [0.008:0.037]   <NA>
    ## 91     [-0.008:0.017]   <NA>
    ## 92     [-0.027:0.016]   <NA>
    ## 93      [0.006:0.033]   <NA>
    ## 94     [-0.008:0.041]   <NA>
    ## 95      [-0.23:0.596]   <NA>
    ## 96     [-0.487:0.399]   <NA>
    ## 97      [0.125:0.433]   <NA>
    ## 98      [0.311:0.572]   <NA>
    ## 99      [0.036:0.296]   <NA>
    ## 100     [0.008:0.319]   <NA>
    ## 101    [-0.005:0.013]   <NA>
    ## 102      [-0.01:0.02]   <NA>
    ## 103     [0.042:0.299]   <NA>
    ## 104      [0.01:0.327]   <NA>
    ## 
    ## $gPartable
    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1    0.17 (0.031)         0.351 (0)     [0.015:0.325]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.247 (0.002)     0.227 (0.003)     [0.089:0.405]
    ## 3   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.266:0.542]
    ## 4     MM_EpiDoc4~MM_EpiDoc2   0.369 (0.016)     0.285 (0.012)     [0.068:0.671]
    ## 5     MM_EpiDoc4~MM_EpiDoc1        0.56 (0)         0.347 (0)     [0.385:0.736]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.102 (0.029)         0.365 (0)      [0.01:0.194]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.494 (0.083)     0.675 (0.028)    [-0.065:1.053]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.077 (0.697)    -0.309 (0.152)    [-0.312:0.466]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.299)     0.086 (0.727)    [-0.205:0.664]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2   0.035 (0.244)     0.001 (0.955)    [-0.024:0.094]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1   0.004 (0.879)     0.028 (0.118)    [-0.052:0.061]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1   0.018 (0.228)      0.009 (0.46)    [-0.011:0.048]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.549 (0.011)     0.922 (0.003)     [0.126:0.973]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.463)      0.01 (0.463)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.046 (0.002)     0.046 (0.002)     [0.017:0.074]
    ## 19           MM_EpiDoc4~Age   0.018 (0.011)         0.031 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age       0 (0.914)     -0.005 (0.29)    [-0.003:0.004]
    ## 21           MM_EpiDoc1~Age       0.032 (0)         0.049 (0)     [0.023:0.041]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.202:0.777]
    ## 24    Dep_EpiDoc1~Education   0.448 (0.002)         0.776 (0)     [0.161:0.735]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education          0 (NA)            0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education          0 (NA)            0 (NA)             [0:0]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4        6.91 (0)          9.37 (0)     [4.942:8.878]
    ## 29 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)     [4.487:7.914]
    ## 30   MM_EpiDoc4~~MM_EpiDoc4       0.971 (0)         0.999 (0)     [0.764:1.178]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.189 (0)         0.401 (0)      [0.087:0.29]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.431 (0)        12.642 (0)    [6.156:10.706]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1       0.701 (0)         0.782 (0)     [0.531:0.871]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1       1.461 (0)         1.293 (0)     [1.042:1.879]
    ## 38            Dep_EpiDoc2~1  -0.261 (0.646)     0.513 (0.401)    [-1.383:0.862]
    ## 39             MM_EpiDoc4~1  -0.362 (0.162)    -0.559 (0.015)    [-0.869:0.145]
    ## 40             MM_EpiDoc2~1   -0.02 (0.756)     0.165 (0.278)    [-0.144:0.105]
    ## 41            Dep_EpiDoc1~1   -0.657 (0.29)    -0.157 (0.819)    [-1.872:0.559]
    ## 42             MM_EpiDoc1~1      -0.785 (0)        -1.335 (0)   [-1.129:-0.441]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.205:0.496]
    ## 2       [0.08:0.374]
    ## 3      [0.208:0.538]
    ## 4      [0.064:0.506]
    ## 5      [0.176:0.519]
    ## 6      [0.177:0.553]
    ## 7      [0.073:1.277]
    ## 8     [-0.731:0.114]
    ## 9     [-0.399:0.572]
    ## 10    [-0.032:0.034]
    ## 11    [-0.007:0.062]
    ## 12    [-0.015:0.034]
    ## 13     [0.315:1.528]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.017:0.074]
    ## 19     [0.017:0.044]
    ## 20    [-0.013:0.004]
    ## 21     [0.038:0.059]
    ## 22             [0:0]
    ## 23       [0.02:0.85]
    ## 24     [0.389:1.162]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27             [0:0]
    ## 28       [7.3:11.44]
    ## 29    [8.104:13.414]
    ## 30     [0.833:1.166]
    ## 31      [0.28:0.521]
    ## 32    [9.241:16.043]
    ## 33      [0.543:1.02]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.809:1.776]
    ## 38    [-0.692:1.719]
    ## 39    [-1.01:-0.107]
    ## 40    [-0.133:0.464]
    ## 41      [-1.5:1.186]
    ## 42   [-1.728:-0.941]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

# GENERATING PLOTs

## prepare for miPLOt

``` r
library(mice)
```

    ## Warning: package 'mice' was built under R version 4.4.3

    ## 
    ## Attaching package: 'mice'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

``` r
data_bal=complete(data.mi,0)
res_dt0<-sem(data=data_bal,
                   model=Model_ED,
                    estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart", 
           sampling.weights='ipw',
         orthogonal=TRUE
  )
```

``` r
lay0=get_layout(
 
   "", "","", 'Dep_EpiDoc1',"","",'',"",  "","Dep_EpiDoc4",
          "", '',"", "","", "", 'Dep_EpiDoc2', "","", "", 
            "",'',"","", "", "","",  "","", "",
    "Age", "","", '',"","",'',"",  "","",
            "", "","", '',"","",'',"",  "","",
    "", "","", '',"","",'',"",  "","",
     "", "","", '',"","",'',"",  "","",
      #     "Education", "","", '',"","",'',"",  "","",
  #          "", "","", '',"","",'',"",  "","",
    "", "","", '',"","",'MM_EpiDoc2',"",  "","",
  "Education",'',"","MM_EpiDoc1", "","",'',"","",  "MM_EpiDoc4",
     #       "", "","", '',"","",'',"",  "","",
       #    "",'',"", '',"","","","", "", '',
           rows = 9)
```

## groups

``` r
dfPar<-get_edges(res_dt0)
p <- prepare_graph(res_dt0, layout = lay0)
dfNPar<-p$nodes
partable<-parameterEstimates.mi(resED,standardized =TRUE)
partable$group<-factor(partable$group,labels=c('special','male','female'))
a=partable%>%mutate(
                      pvalue_temp=ifelse(pvalue<0.001,'<0.001',sprintf('=%.3f',pvalue)),
                      label2=ifelse(is.na(pvalue_temp),'',
                        sprintf('%.3f (p%s)',std.all,pvalue_temp)))%>%
  group_by(rhs,lhs,op)%>%summarise(label=paste(paste(group,': ',label2),collapse=';\n'))%>%
  select(rhs,lhs,op,label)%>%mutate(
    label=ifelse(op=='~~',
                                    ifelse(rhs==lhs,'',label),
                                    ifelse(op=='~1','',label))
  )%>%mutate(
    label=ifelse(label=='male :  ;\nfemale :  ','',label)
  )
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by rhs, lhs, and op.
    ## ℹ Output is grouped by rhs and lhs.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(rhs, lhs, op))` for per-operation grouping
    ##   (`?dplyr::dplyr_by`) instead.

``` r
dfPar2<-left_join(dfPar%>%select(rhs,lhs,op,arrow,est_sig_std),
                    a)
```

    ## Joining with `by = join_by(rhs, lhs, op)`

``` r
b<-prepare_graph(res_dt0, layout = lay0,
                 rect_width=3,rect_height=3,
                 spacing_y=3)%>%
  edit_graph( {label = dfPar2$label},element='edges') %>%
  edit_graph({label_location = 0.3})%>%
  edit_graph({label_size = 3},element='edges')%>%
  color_var(NA) 

c <- b$edges[b$edges$label!='', ]
c <- c[c$pval_std<0.10, ]
c=c%>%mutate(color=ifelse(pval_std>0.05,'grey',ifelse(est_std>0,'blue','red')),
             label_color=ifelse(pval_std>0.05,'grey',
                                ifelse(est_std>0,'blue','red')),
             label_location=ifelse(from=='Age',runif(n(),min=0.1,0.2),
                                   ifelse(from=="Education",runif(n(),min=0.1,0.4),0.5)))
b$edges<-c
b <- edit_graph(b, 
                element = "edges", 
                label_alpha = 0,
                label_fill = "transparent")
b$nodes$fill='steelblue1'
b$nodes$color='dodgerblue3'
b$nodes$label_fill='transparent'
b%>%
  plot(size = c(15, 10))
```

    ## Some edges involve nodes with argument 'show = FALSE'. These were dropped.

![](5_plot_git_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave(filename=paste0('./Results/finalModel/laggedSEM_MGA.png'), width = 35, height = 15, units = "cm")
```

## full model

``` r
dfPar<-get_edges(res_dt0)

b<-prepare_graph(res_dt0, layout = lay0,
                 rect_width=3,rect_height=3,
                 spacing_y=3,
                 hide_mean = TRUE)%>%
  edit_graph({label_location = 0.3})%>%
  edit_graph({label_size = 3},element='edges')%>%
  color_var(NA) 
c <- b$edges[b$edges$label!='', ]
c[c$lhs==c$rhs, 'label']<-''
c[c$lhs==c$rhs, 'show']<-FALSE
c <- c[c$pval_std<0.10, ]
c=c%>%mutate(color=ifelse(pval_std>0.05,'grey',ifelse(est_std>0,'blue','red')),
             label_color=ifelse(pval_std>0.05,'grey',
                                ifelse(est_std>0,'blue','red')),
             label_location=ifelse(from=='Age',runif(n(),min=0.1,0.2),
                                   ifelse(from=="Education",runif(n(),min=0.1,0.4),0.5)),
             label=ifelse(pval_std<0.001, paste0(est_std,'(p=<0.001)'),
               paste0(est_std,'(p=',pval_std,')')))
b$edges<-c
b <- edit_graph(b, 
                element = "edges", 
                label_alpha = 0,
                label_fill = "transparent")
b$nodes$fill='steelblue1'
b$nodes$color='dodgerblue3'
b$nodes$label_fill='transparent'
b%>%
  plot(size = c(15, 10))
```

    ## Some edges involve nodes with argument 'show = FALSE'. These were dropped.

![](5_plot_git_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave(filename=paste0('./Results/finalModel/laggedSEM_avgModel.png'), width = 35, height = 15, units = "cm")
```

``` r
devtools::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", :
    ## running command '"quarto"
    ## TMPDIR=C:/Users/tomoe/AppData/Local/Temp/RtmpIjeTvi/file43883b53542c -V' had
    ## status 1

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.4.2 (2024-10-31 ucrt)
    ##  os       Windows 11 x64 (build 26200)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United Kingdom.utf8
    ##  ctype    English_United Kingdom.utf8
    ##  tz       America/Sao_Paulo
    ##  date     2026-06-01
    ##  pandoc   3.6.3 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   NA @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package         * version  date (UTC) lib source
    ##  abind             1.4-8    2024-09-12 [1] CRAN (R 4.4.1)
    ##  backports         1.5.0    2024-05-23 [1] CRAN (R 4.4.0)
    ##  boot              1.3-31   2024-08-28 [1] CRAN (R 4.4.2)
    ##  broom             1.0.9    2025-07-28 [1] CRAN (R 4.4.3)
    ##  cachem            1.1.0    2024-05-16 [1] CRAN (R 4.4.3)
    ##  car               3.1-3    2024-09-27 [1] CRAN (R 4.4.3)
    ##  carData           3.0-5    2022-01-06 [1] CRAN (R 4.4.3)
    ##  checkmate         2.3.3    2025-08-18 [1] CRAN (R 4.4.3)
    ##  cli               3.6.5    2025-04-23 [1] CRAN (R 4.4.3)
    ##  coda              0.19-4.1 2024-01-31 [1] CRAN (R 4.4.3)
    ##  codetools         0.2-20   2024-03-31 [1] CRAN (R 4.4.2)
    ##  CompQuadForm      1.4.4    2025-07-13 [1] CRAN (R 4.4.3)
    ##  data.table        1.17.8   2025-07-10 [1] CRAN (R 4.4.3)
    ##  dbscan            1.2.3    2025-08-20 [1] CRAN (R 4.4.3)
    ##  devtools          2.5.0    2026-03-14 [1] CRAN (R 4.4.3)
    ##  digest            0.6.37   2024-08-19 [1] CRAN (R 4.4.3)
    ##  dplyr           * 1.2.1    2026-04-03 [1] CRAN (R 4.4.2)
    ##  ellipsis          0.3.2    2021-04-29 [1] CRAN (R 4.4.3)
    ##  evaluate          1.0.4    2025-06-18 [1] CRAN (R 4.4.3)
    ##  farver            2.1.2    2024-05-13 [1] CRAN (R 4.4.3)
    ##  fastDummies       1.7.5    2025-01-20 [1] CRAN (R 4.4.3)
    ##  fastmap           1.2.0    2024-05-15 [1] CRAN (R 4.4.3)
    ##  foreach           1.5.2    2022-02-02 [1] CRAN (R 4.4.3)
    ##  Formula           1.2-5    2023-02-24 [1] CRAN (R 4.4.0)
    ##  fs                2.0.1    2026-03-24 [1] CRAN (R 4.4.3)
    ##  future            1.70.0   2026-03-14 [1] RSPM
    ##  future.apply      1.20.0   2025-06-06 [1] CRAN (R 4.4.3)
    ##  generics          0.1.4    2025-05-09 [1] CRAN (R 4.4.3)
    ##  ggplot2         * 4.0.2    2026-02-03 [1] CRAN (R 4.4.3)
    ##  glmnet            4.1-10   2025-07-17 [1] CRAN (R 4.4.3)
    ##  globals           0.19.1   2026-03-13 [1] RSPM
    ##  glue              1.8.0    2024-09-30 [1] CRAN (R 4.4.3)
    ##  gsubfn            0.7      2018-03-16 [1] CRAN (R 4.4.3)
    ##  gtable            0.3.6    2024-10-25 [1] CRAN (R 4.4.3)
    ##  htmltools         0.5.8.1  2024-04-04 [1] CRAN (R 4.4.3)
    ##  httr              1.4.7    2023-08-15 [1] CRAN (R 4.4.3)
    ##  igraph            2.1.4    2025-01-23 [1] CRAN (R 4.4.3)
    ##  iterators         1.0.14   2022-02-05 [1] CRAN (R 4.4.3)
    ##  jomo              2.7-6    2023-04-15 [1] CRAN (R 4.4.3)
    ##  knitr             1.50     2025-03-16 [1] CRAN (R 4.4.3)
    ##  labeling          0.4.3    2023-08-29 [1] CRAN (R 4.4.0)
    ##  lattice           0.22-6   2024-03-20 [1] CRAN (R 4.4.2)
    ##  lavaan          * 0.6-19   2024-09-26 [1] CRAN (R 4.4.3)
    ##  lavaan.mi       * 0.1-0    2025-03-10 [1] CRAN (R 4.4.3)
    ##  lifecycle         1.0.5    2026-01-08 [1] CRAN (R 4.4.3)
    ##  listenv           0.9.1    2024-01-29 [1] CRAN (R 4.4.3)
    ##  lme4              1.1-37   2025-03-26 [1] CRAN (R 4.4.3)
    ##  magrittr          2.0.5    2026-04-04 [1] RSPM
    ##  MASS              7.3-61   2024-06-13 [1] CRAN (R 4.4.2)
    ##  Matrix            1.7-1    2024-10-18 [1] CRAN (R 4.4.2)
    ##  memoise           2.0.1    2021-11-26 [1] CRAN (R 4.4.3)
    ##  mice            * 3.18.0   2025-05-27 [1] CRAN (R 4.4.3)
    ##  minqa             1.2.8    2024-08-17 [1] CRAN (R 4.4.3)
    ##  mitml             0.4-5    2023-03-08 [1] CRAN (R 4.4.3)
    ##  mnormt            2.1.1    2022-09-26 [1] CRAN (R 4.4.0)
    ##  MplusAutomation   1.2      2025-09-02 [1] CRAN (R 4.4.3)
    ##  mvtnorm           1.3-3    2025-01-10 [1] CRAN (R 4.4.3)
    ##  nlme              3.1-166  2024-08-14 [1] CRAN (R 4.4.2)
    ##  nloptr            2.2.1    2025-03-17 [1] CRAN (R 4.4.3)
    ##  nnet              7.3-19   2023-05-03 [1] CRAN (R 4.4.2)
    ##  nonnest2          0.5-8    2024-08-28 [1] CRAN (R 4.4.3)
    ##  pan               1.9      2023-08-21 [1] CRAN (R 4.4.3)
    ##  pander            0.6.6    2025-03-01 [1] CRAN (R 4.4.3)
    ##  parallelly        1.45.1   2025-07-24 [1] CRAN (R 4.4.3)
    ##  pbivnorm          0.6.0    2015-01-23 [1] CRAN (R 4.4.0)
    ##  pillar            1.11.1   2025-09-17 [1] CRAN (R 4.4.3)
    ##  pkgbuild          1.4.8    2025-05-26 [1] CRAN (R 4.4.3)
    ##  pkgconfig         2.0.3    2019-09-22 [1] CRAN (R 4.4.3)
    ##  pkgload           1.5.1    2026-04-01 [1] CRAN (R 4.4.2)
    ##  plyr              1.8.9    2023-10-02 [1] CRAN (R 4.4.3)
    ##  progressr         0.15.1   2024-11-22 [1] CRAN (R 4.4.3)
    ##  proto             1.0.0    2016-10-29 [1] CRAN (R 4.4.3)
    ##  psych             2.5.6    2025-06-23 [1] CRAN (R 4.4.3)
    ##  purrr             1.2.1    2026-01-09 [1] CRAN (R 4.4.3)
    ##  quadprog          1.5-8    2019-11-20 [1] CRAN (R 4.4.0)
    ##  R6                2.6.1    2025-02-15 [1] CRAN (R 4.4.3)
    ##  ragg              1.5.2    2026-03-23 [1] CRAN (R 4.4.3)
    ##  RANN              2.6.2    2024-08-25 [1] CRAN (R 4.4.3)
    ##  rbibutils         2.3      2024-10-04 [1] CRAN (R 4.4.3)
    ##  RColorBrewer      1.1-3    2022-04-03 [1] CRAN (R 4.4.0)
    ##  Rcpp              1.1.0    2025-07-02 [1] CRAN (R 4.4.3)
    ##  Rdpack            2.6.4    2025-04-09 [1] CRAN (R 4.4.3)
    ##  reformulas        0.4.1    2025-04-30 [1] CRAN (R 4.4.3)
    ##  renv              1.1.5    2025-07-24 [1] CRAN (R 4.4.3)
    ##  rlang             1.2.0    2026-04-06 [1] CRAN (R 4.4.2)
    ##  rmarkdown         2.29     2024-11-04 [1] CRAN (R 4.4.3)
    ##  rpart             4.1.23   2023-12-05 [1] CRAN (R 4.4.2)
    ##  rstudioapi        0.17.1   2024-10-22 [1] CRAN (R 4.4.3)
    ##  S7                0.2.1    2025-11-14 [1] CRAN (R 4.4.3)
    ##  sandwich          3.1-1    2024-09-15 [1] CRAN (R 4.4.3)
    ##  scales            1.4.0    2025-04-24 [1] CRAN (R 4.4.3)
    ##  sessioninfo       1.2.3    2025-02-05 [1] CRAN (R 4.4.3)
    ##  shape             1.4.6.1  2024-02-23 [1] CRAN (R 4.4.0)
    ##  survival          3.7-0    2024-06-05 [1] CRAN (R 4.4.2)
    ##  systemfonts       1.3.1    2025-10-01 [1] CRAN (R 4.4.3)
    ##  texreg            1.39.4   2024-07-24 [1] CRAN (R 4.4.3)
    ##  textshaping       1.0.4    2025-10-10 [1] CRAN (R 4.4.3)
    ##  tibble            3.3.1    2026-01-11 [1] CRAN (R 4.4.3)
    ##  tidyr             1.3.2    2025-12-19 [1] CRAN (R 4.4.3)
    ##  tidyselect        1.2.1    2024-03-11 [1] CRAN (R 4.4.3)
    ##  tidySEM         * 0.2.9    2025-07-30 [1] CRAN (R 4.4.3)
    ##  usethis           3.2.1    2025-09-06 [1] CRAN (R 4.4.3)
    ##  vctrs             0.7.2    2026-03-21 [1] CRAN (R 4.4.3)
    ##  withr             3.0.2    2024-10-28 [1] CRAN (R 4.4.3)
    ##  xfun              0.53     2025-08-19 [1] CRAN (R 4.4.3)
    ##  xtable            1.8-4    2019-04-21 [1] CRAN (R 4.4.3)
    ##  yaml              2.3.10   2024-07-26 [1] CRAN (R 4.4.3)
    ##  zoo               1.8-14   2025-04-10 [1] CRAN (R 4.4.3)
    ## 
    ##  [1] D:/1_Data/RutePortugal/renv/library/windows/R-4.4/x86_64-w64-mingw32
    ##  [2] C:/Users/tomoe/AppData/Local/R/cache/R/renv/sandbox/windows/R-4.4/x86_64-w64-mingw32/6698a5f3
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
