# Incomplete Beta Function

This is [CodePlea's implementation of the regularized incomplete beta
function in C](https://codeplea.com/incomplete-beta-function-c),
modified to make it easy to call from Dyalog APL.


## Building the shared library

```
$ cmake . -DCMAKE_BUILD_TYPE=Release
$ make
```


## Using it from APL

```APL
      lib←'libincbeta.so' ⍝ varies according to platform
      ⎕NA lib,'|incbeta F8 F8 =F8[] P'                                
      A←1 ⋄ B←3 ⋄ 4⍕X←?10/0
 0.9403 0.9214 0.4226 0.2387 0.9489 0.3767 0.6916 0.8186 0.5756 0.2557
      4⍕incbeta A B X (≢X)
 0.9998 0.9995 0.8075 0.5588 0.9999 0.7578 0.9707 0.9940 0.9236 0.5876
```
