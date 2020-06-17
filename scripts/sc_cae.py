#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 21:41:04 2019

@author: drmz
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

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
x = tf.keras.Input(shape=(dim_orig,))
h = tf.keras.layers.Dense(512, activation='relu', use_bias=True)(x)
h = tf.keras.layers.Dense(128, activation='relu', use_bias=True)(h)
h = tf.keras.layers.Dense(32, activation='relu', use_bias=True)(h)

z_mean = tf.keras.layers.Dense(dim_lat, activation=None)(h)
z_log_var = tf.keras.layers.Dense(dim_lat, activation=None)(h)

# Reparameterize
z = tf.keras.layers.Lambda(sampling, output_shape=(dim_lat,))([z_mean, z_log_var])

# Decoder
h_decoded = tf.keras.layers.Dense(32, activation='relu')(z)
h_decoded = tf.keras.layers.Dense(128, activation='relu')(h_decoded)
h_decoded = tf.keras.layers.Dense(512, activation='relu')(h_decoded)
x_decoded_mean = tf.keras.layers.Dense(dim_orig,activation='linear')(h_decoded)


vae = tf.keras.Model(x, x_decoded_mean)

# Encoding
enc = tf.keras.Model(x, z_mean)

# Loss
def vae_loss(x, x_decoded_mean):
    mse_loss = dim_orig*tf.keras.losses.mean_squared_error(x, x_decoded_mean)
    kl_loss = 1 + z_log_var - tf.keras.backend.square(z_mean) - tf.keras.backend.exp(z_log_var)
    kl_loss = tf.keras.backend.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    
    local_loss = tf.keras.backend.mean(mse_loss + kl_loss)
    return local_loss

# Optimizer
opt_adam = tf.keras.optimizers.Adam(lr=0.0001)
vae.compile(optimizer=opt_adam, loss=vae_loss)
model_history = vae.fit(expr_data, expr_data, epochs=50, batch_size=batch_size, validation_split=0.2, shuffle=True)

# Save Embedding
c_encoded = enc.predict(expr_data)

np.savetxt(out_file, c_encoded)
tf.keras.backend.clear_session()
