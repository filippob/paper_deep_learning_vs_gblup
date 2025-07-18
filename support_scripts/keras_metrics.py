"""
Created on Thu Jan 20 17:14:30 2022

@author: Nelson
"""

""" A collection of custom metrics for keras """

import numpy as np
import tensorflow as tf
import keras.backend as KB

#pearson's correlation
def pearson(x, y):
  #formula https://www.statology.org/pearson-correlation-coefficient/

  #NUMERATOR: SUM[ (x - x_mean) * (y - y_mean)]
  #DENOMINATOR: SQRT(A) * SQRT(B)
  #with:
  #  A = SUM[ (x - x_mean) ^ 2]
  #  B = SUM[ (y - y_mean) ^ 2]

  num = tf.math.reduce_sum( (x - tf.math.reduce_mean(x)) * (y - tf.math.reduce_mean(y)) )
  den = tf.math.sqrt(tf.math.reduce_sum((x - tf.math.reduce_mean(x)) ** 2)) * tf.math.sqrt(tf.math.reduce_sum((y - tf.math.reduce_mean(y)) ** 2))

  return(num / den)

#Root Mean Square Error, to ease comparison with GROAN
def rmse(x, y):
  return tf.math.reduce_mean(tf.math.sqrt((x - y) ** 2))

## NDCG: normalised discounted cumulative gain
## 1) basic version to work with arrays (numpy))
def ndcg(y, y_hat, k):
    
    n = len(y)
    ## select the k top examples
    nk = np.round(k*n).astype(int)

    ## decreasing order: use slicing
    ## arr[start:end:step]
    y_inds = np.argsort(y)[::-1] ##revert argsort to get decreasing order
    y_sort_y = y[y_inds]
    y_hat_inds = np.argsort(y_hat)[::-1]
    y_sort_y_hat = y[y_hat_inds]
    
    seq = np.arange(1.0, nk+1)
    d = 1/np.log2(seq+1)
    
    ## subset arrays
    sliced_y_hat = y_sort_y_hat[0:nk]
    sliced_y = y_sort_y[0:nk]
    
    num = sum(sliced_y_hat*d)
    den = sum(sliced_y*d)

    temp = num/den
    
    return(temp)


## 2) version for tensors (Keras) [in progress]
def ndcg_tf(y, y_hat, k):
    
    #n = len(y)
    kt = tf.convert_to_tensor(k, dtype=tf.float32)
    nt = tf.convert_to_tensor(len(y), dtype=tf.int32) # number of predicted examples
    nt = tf.cast(nt, dtype=tf.float32)
    ## select the k top examples
    nk = tf.math.round(kt*nt)
    #nk_int = round(n*k)
    nk_int = tf.cast(nk, dtype=tf.int32)
    
    y_inds = tf.argsort(y, direction='DESCENDING')
    y_sort_y = tf.gather(y, y_inds)
    y_hat_inds = tf.argsort(y_hat, direction='DESCENDING')
    y_sort_y_hat = tf.gather(y, y_hat_inds)
    
    seq = tf.range(1.0, nk+1)
    d = KB.log(seq+1)
    k2 = KB.log(2.0)
    d = k2/d ## 1/(d/k2) --> 1* k2/d
    
    print('shape of tensor y', tf.shape(y))
    y_print = KB.eval(y)
    print(y_print)
    print('shape of tensor y_inds', tf.shape(y_inds))
    print(y_inds)
    
    ## numpy array syntax: does not work with Keras eager execution
    ## num = KB.sum(y_sort_y_hat[0:nk_int]*d)
    ## den = KB.sum(y_sort_y[0:nk_int]*d)
    
    print('shape of sorted tensor num', tf.shape(y_sort_y_hat))
    print('shape of sorted tensor num', y_sort_y_hat.get_shape())
    print(y_sort_y_hat)
    print('shape of sorted tensor den', tf.shape(y_sort_y))
    print('shape of sorted tensor den', y_sort_y.get_shape())
    print(y_sort_y)
    
    ## tensor flow slice syntax: for Keras with eager execution
    sliced_y_hat = tf.slice(y_sort_y_hat, [0], [nk_int])
    sliced_y = tf.slice(y_sort_y, [0], [nk_int])
    
    num = KB.sum(sliced_y_hat*d)
    den = KB.sum(sliced_y*d)

    temp = num/den
    
    return(temp)

def ndcg_25(y, yhat):
    
    return(ndcg_tf(y, yhat, 0.25))
    

def ndcg_50(y, yhat):
    
    return(ndcg_tf(y, yhat, 0.50))

def ndcg_100(y, yhat):
    
    return(ndcg_tf(y, yhat, 1.0))
    
def ndcg_nok(y, y_hat):
    
    k = 1.0
    kt = tf.convert_to_tensor(k, dtype=tf.float32)
    nt = tf.convert_to_tensor(len(y), dtype=tf.int32) # number of predicted examples
    nt = tf.cast(nt, dtype=tf.float32)
    ## select the k top examples
    nk = tf.math.round(kt*nt)
    
    y_inds = tf.argsort(y, direction='DESCENDING')
    y_sort_y = tf.gather(y, y_inds)
    y_hat_inds = tf.argsort(y_hat, direction='DESCENDING')
    y_sort_y_hat = tf.gather(y, y_hat_inds)
    
    seq = tf.range(1.0, nk+1)
    d = KB.log(seq+1)
    k2 = KB.log(2.0)
    d = k2/d ## 1/(d/k2) --> 1* k2/d
    
    num = KB.sum(y_sort_y_hat*d)
    den = KB.sum(y_sort_y*d)

    temp = num/den
    
    return(temp)
