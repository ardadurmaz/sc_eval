#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 21:41:04 2019

@author: drmz
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import tensorflow as tf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt


data_file = sys.argv[1]
out_file = sys.argv[2]

# Load Data
expr_data = np.loadtxt(data_file, dtype = 'float64')

# VAE
batch_size = 100

dim_orig = expr_data.shape[1]
dim_lat = 3


def sampling(args):
    z_mean, z_log_var = args
    local_batch = tf.keras.backend.shape(z_mean)[0]
    local_dim = tf.keras.backend.int_shape(z_mean)[1]
    eps = tf.keras.backend.random_normal(shape=(local_batch, local_dim))
    
    return z_mean + tf.keras.backend.exp(0.5 * z_log_var)*eps

# Encoder
x = keras.Input(shape = (dim_orig,))
h = keras.layers.Dense(512, activation='relu')(x)
h = keras.layers.Dense(256, activation='relu')(h)
h = keras.layers.Dense(128, activation='relu')(h)
h = keras.layers.Dense(64, activation='relu')(h)

z_mean = keras.layers.Dense(dim_lat)(h)
z_log_var = keras.layers.Dense(dim_lat)(h)

# Reparameterize
z = keras.layers.Lambda(sampling, output_shape=(dim_lat,))([z_mean, z_log_var])

# Decoder
h_decoded = keras.layers.Dense(64, activation='relu')(z)
h_decoded = keras.layers.Dense(128, activation='relu')(h_decoded)
h_decoded = keras.layers.Dense(256, activation='relu')(h_decoded)
h_decoded = keras.layers.Dense(512, activation='relu')(h_decoded)
x_decoded_mean = keras.layers.Dense(dim_orig,activation='relu')(h_decoded)


vae = keras.Model(x, x_decoded_mean)
encoder = keras.Model(x, z_mean)

# Loss
def vae_loss(x, x_decoded_mean):
    #ent_loss = dim_orig*tf.keras.losses.binary_crossentropy(x, x_decoded_mean)
    mse_loss = dim_orig*tf.keras.losses.mean_squared_error(x, x_decoded_mean)
    kl_loss = 1 + z_log_var - tf.keras.backend.square(z_mean) - tf.keras.backend.exp(z_log_var)
    kl_loss = tf.keras.backend.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    
    local_loss = tf.keras.backend.mean(mse_loss + kl_loss)
    return local_loss

# Optimizer
#rms_prop = keras.optimizers.RMSprop(lr=0.00001)
opt_adam = keras.optimizers.Adam(lr=0.0001)
vae.compile(optimizer=opt_adam, loss=vae_loss)
vae.fit(expr_data, expr_data, epochs=100, batch_size=batch_size)


## Save Embedding  ##
c_encoded = encoder.predict(expr_data)
plt.scatter(c_encoded[::,0], c_encoded[::,2])

np.savetxt(out_file, c_encoded)
keras.backend.clear_session()
