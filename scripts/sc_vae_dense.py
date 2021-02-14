#!/usr/bin/env python3

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import tensorflow as tf

tf.compat.v1.disable_eager_execution()


def sampling(args):
    z_mean, z_log_var = args
    local_batch = tf.keras.backend.shape(z_mean)[0]
    local_dim = tf.keras.backend.shape(z_mean)[1]
    eps = tf.keras.backend.random_normal(shape=(local_batch, local_dim))
    return z_mean + tf.keras.backend.exp(0.5 * z_log_var)*eps

def vae_loss_mse(args):
    z_mean, z_log_var = args
    kl_loss = 1.0 + z_log_var - tf.keras.backend.exp(z_log_var) - tf.keras.backend.square(z_mean) 
    kl_loss = tf.keras.backend.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    
    def loss(y_true, y_pred):
        mse_loss = tf.keras.backend.sum(tf.keras.backend.square(y_pred-y_true), axis=-1)
        return tf.keras.backend.mean(mse_loss + kl_loss)
    return loss

def scvae(data):
    BATCH_SIZE = 100
    DIM_ORIG = data.shape[1]
    LAT_DIM = 2
    NUM_EPOCHS = 300
    
    
    encoder_input = tf.keras.layers.Input(shape=(DIM_ORIG,))
    h = tf.keras.layers.Dense(1024, activation='relu')(encoder_input)
    h = tf.keras.layers.Dense(512, activation='relu')(h)
    h = tf.keras.layers.Dense(256, activation='relu')(h)
    
    z_mean = tf.keras.layers.Dense(LAT_DIM)(h)
    z_log_var = tf.keras.layers.Dense(LAT_DIM)(h)
    z = tf.keras.layers.Lambda(sampling, output_shape=(LAT_DIM,), name='z')([z_mean, z_log_var])
    
    h_decoded = tf.keras.layers.Dense(256, activation='relu')(z)
    h_decoded = tf.keras.layers.Dense(512, activation='relu')(h_decoded)
    h_decoded = tf.keras.layers.Dense(1024, activation='relu')(h_decoded)
    encoder_output = tf.keras.layers.Dense(DIM_ORIG)(h_decoded)
         
    vae = tf.keras.Model(inputs=encoder_input, outputs=encoder_output)
    encoder_model = tf.keras.Model(encoder_input, z_mean, name='encoder')
    
    # Callback
    ear_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, mode='auto')
    
    # Fit
    vae.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate=0.0001), loss=vae_loss_mse([z_mean, z_log_var]))
    vae.summary()
    
    vae.fit(data, data, epochs=NUM_EPOCHS, batch_size=BATCH_SIZE, shuffle=True, callbacks=[ear_stop], validation_split=0.2)
    
    # Embedding
    embedd=encoder_model.predict(data)

    tf.keras.backend.clear_session()
    return embedd
    
if __name__ == '__main__':
    # Read Data
    local_data = np.loadtxt('Expression.tsv', dtype='float64')
    
    # Run VAE
    embedd = scvae(local_data)
    
    # Save
    np.savetxt(fname='VAE_Embedding.tsv', X=embedd, delimiter='\t')
