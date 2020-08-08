# -*- coding: utf-8 -*-
"""
Created on %25 Apr 2018(12:05pm)

@author: %A.Nishanth C00294860
"""
import tensorflow as tf
import os
import numpy as np
import pickle
#%% then load the training samples subgroups from the directory
"""
Here use every think as training sample inorder to identify the MOTIFs among them

We can do in pleanty of ways

    Just choose the PDB_ids subgroup as sets
    use 70 % of the PDB_ids_groups as training
    use 20 % of the PDB_ids_groups as testing
    use 10 % of the PDB_ids_groups as validation
    
1-Way:
    
    Randomize all groups with label and 
    
    Train with the training set the
    
    >>  loss is normal subgroup missclasification

2-Way:
    
    PDBgroup set is whole one set not-seperable
    
    Train with the training set 
    Then the subgroup is classified as ONGO or TSG MOTIF  
        then the over all class(ONGO or TSG) is choosen from the mjority vote of the subgroups
        
    >>  loss is from the PDB_class speicifiaction error
    
>>> ---------------------------------------------------------------------------

The 2-way has some issues
    Since Each PDB has atleat 100 subclasses
    Training is expensive for computation
"""
def accuracy(predictions, batch_labels):
    acc = 0
    for i in range(0,len(batch_labels)):
        if batch_labels[i]==predictions[i]:
            acc = acc + 1
    return 100*acc/len(batch_labels)
#%% load the way_1 training set
working_dir = "saving directory.../Deep_learning_data_initialization"
os.chdir('/')
os.chdir(working_dir)
train = np.load('shuffled_all_data.npy') 
train_label = np.load('shuffled_all_label.npy') 
#%% neural netwoek formation starts from here
"""
Neural Network design 

    Here 64 is choosen inorder to extract the 
    (atleast 4 amino acids needed to form a substructure)
                4 aminoacids with 16 properties  
                4 x 16 = 64
                
    weight array1: 1715 x 64
                followed by relu inorder to make nonlinearity
    weight array2: 
                64 x 1 to findout ONGO or TSG

"""

batch_size = 500 # since one PDB structure contain atleast 100 sub structures 
num_steps = 4000 # almost gothrough 5 times each training sample
learning_rate = 1e-2


graph = tf.Graph()

with graph.as_default():
  # Input data.
  tf_train_dataset = tf.placeholder(tf.float32, shape=(batch_size,1715),name="tf_train_dataset")
  tf_train_labels = tf.placeholder(tf.float32, shape=(batch_size, 1),name="tf_train_labels")

  layer1_weights = tf.Variable(tf.truncated_normal([1715,64], stddev=0.1))
  layer1_biases = tf.Variable(tf.zeros([64]))
  layer2_weights = tf.Variable(tf.truncated_normal([64,1], stddev=0.1))
  layer2_biases = tf.Variable(tf.constant(1.0, shape=[1]))
 
  # for dropout
#  keep_prob = tf.placeholder(tf.float32,name="keep_prob")
 
  # Model.
  def model(data):
    #layer 1
    layer_1 = tf.nn.relu(tf.matmul(data, layer1_weights) + layer1_biases, name='hidden_1')
    return tf.matmul(layer_1,  layer2_weights) + layer2_biases

  # Training computation.
  logits = model(tf_train_dataset)
#  loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(labels=tf_train_labels, logits=logits))
  loss = tf.losses.log_loss(tf_train_labels,tf.nn.softmax(logits))
  # Optimizer.
#  optimizer = tf.train.GradientDescentOptimizer(learning_rate).minimize(loss)
  optimizer = tf.train.AdamOptimizer(learning_rate).minimize(loss)
  # Predictions for the training, validation, and test data.
  train_prediction = tf.nn.softmax(logits,name="train_prediction")
  # Add ops to save and restore all the variables.
  saver = tf.train.Saver()
#%% 
accuracy_vector_train=[]#to save the train accuracy details
loss_vector=[]
final_batch_index = []
final_batch_corresponding_prediction= []
with tf.Session(graph=graph) as session:
    tf.global_variables_initializer().run()

    print('')
    print(" ")    
    print('Initialized')   
    print('')
    print('')
    for step in range(num_steps):
        offset = (step * batch_size) % (train_label.shape[0] - batch_size)
#        offset=0
        batch_data = train[offset:(offset + batch_size), :]
        batch_labels = train_label[offset:(offset + batch_size)]
        feed_dict = {tf_train_dataset : batch_data, tf_train_labels : batch_labels}
#        feed_dict = {tf_train_dataset : batch_data, tf_train_labels : batch_labels}
        _, l, predictions = session.run([optimizer, loss, train_prediction], feed_dict=feed_dict)
        if (step % 50 == 0):
            print('')
            print('step---: ',step)
            print('Minibatch accuracy: %.1f%%' % accuracy(predictions, batch_labels))
            print('Minibatch loss    : ', l*(10^6))
            accuracy_vector_train.append(accuracy(predictions, batch_labels)) #to save the accuracy details
            loss_vector.append(l)
    
    for step in range(0,765):                
        offset = (step * batch_size) % (train_label.shape[0] - batch_size)
        batch_data = train[offset:(offset + batch_size), :]
        batch_labels = train_label[offset:(offset + batch_size), :]
        feed_dict = {tf_train_dataset : batch_data, tf_train_labels : batch_labels}
    #        feed_dict = {tf_train_dataset : batch_data, tf_train_labels : batch_labels}
        predictions = session.run([train_prediction], feed_dict=feed_dict)
        final_batch_index.append(range(offset,(offset + batch_size)))
#        print("run chk------")
        final_batch_corresponding_prediction.append(predictions)
#    print("save check------")
#to access        
#        final_batch_corresponding_prediction[0][0][499][0]
    os.chdir('/')
    os.chdir(working_dir)
    save = saver.save(session, ''.join([working_dir,"/","train_model_1",".ckpt"]),global_step=1)  
    pickle.dump(final_batch_index, open("final_batch_index.p", "wb")) 
    pickle.dump(final_batch_corresponding_prediction, open("final_batch_corresponding_prediction.p", "wb")) 
    pickle.dump(accuracy_vector_train, open("accuracy_vector_train.p", "wb")) 
    pickle.dump(loss_vector, open("loss_vector.p", "wb")) 
#%% load the results and check the performance
final_batch_index = pickle.load(open("final_batch_index.p", "rb")) 
final_batch_corresponding_prediction = pickle.load(open("final_batch_corresponding_prediction.p", "rb"))  
accuracy_vector_train = pickle.load(open("accuracy_vector_train.p", "rb"))  
loss_vector = pickle.load(open("loss_vector.p", "rb"))  
