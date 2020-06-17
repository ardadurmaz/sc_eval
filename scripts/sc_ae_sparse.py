#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:05:14 2019

@author: drmz
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import tensorflow as tf


class SparseL(tf.keras.layers.Layer):
    def __init__(self, units=None, reg_mat=None, dtype=tf.float32):
        super(SparseL, self).__init__()
        self.units = units
        self.group_w = reg_mat
    
    def build(self, input_shape):
        self.kernel = self.add_weight('kernel', shape=[int(input_shape[-1]), self.units], initializer='glorot_uniform', dtype=self.dtype, trainable=True)
        self.bias = self.add_weight('bias', shape=[self.units,], initializer='glorot_uniform', dtype=self.dtype, trainable=True)
    
    def call(self, inputs):
        return tf.math.sigmoid(tf.matmul(inputs, tf.math.multiply(self.group_w, self.kernel)) + self.bias)


def scae(in_mat=None, out_mat=None, w1=None, w2=None):
    BATCH_SIZE = 500
    DIM_IN = in_mat.shape[1]
    DIM_OUT = out_mat.shape[1]
    DIM_LAT = w1.shape[1]
    NUM_EPOCHS = 50
        
    encoder_input = tf.keras.layers.Input(shape=(DIM_IN,))
    h = SparseL(50, reg_mat=w1)(encoder_input)
    encoder_output = SparseL(DIM_OUT, reg_mat=w2)(h)
    
    encoder_model = tf.keras.Model(encoder_input, h, name='encoder')
    scae_model = tf.keras.Model(inputs=encoder_input, outputs=encoder_output)
    
    scae_model.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate=0.0001), loss='mse')
    scae_model.summary()
    
    ## Fit ##
    scae_model.fit(in_mat, out_mat, epochs=NUM_EPOCHS, batch_size=BATCH_SIZE)
    
    
    ## Encoding ##
    embedd = encoder_model.predict(in_mat)
    tf.keras.backend.clear_session()
    return(embedd)
  
if __name__ == '__main__':
    # Read Data
    in_mat = np.loadtxt('/home/drmz/Research/sc_eval/InputMat.tsv', dtype=np.float32)
    out_mat = np.loadtxt('/home/drmz/Research/sc_eval/OutputMat.tsv', dtype=np.float32)
    w1 = np.loadtxt('/home/drmz/Research/sc_eval/W1.tsv', dtype=np.float32)
    w2 = np.loadtxt('/home/drmz/Research/sc_eval/W2.tsv', dtype=np.float32)
    
    # Run VAE
    embedd = scae(local_data..)
    print("Saving Embeddings..")
    np.savetxt(fname='/home/drmz/Research/sc_eval/TestEmbedding.tsv', X=embedd, delimiter='\t')    

