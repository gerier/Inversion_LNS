from copy import deepcopy

def model_actualisation(model, kernel, step):
    # operation of gradient descent
    new_model = model - step * kernel 
    return new_model 

def get_derivative(model):
    # get the derivative of the wanted function (cost function)
    kernel = 2 * model
    return kernel

def inversion(model_0, step=0.1):

    # initialisation
    model = deepcopy(model_0)
    it = 0
    criteria = it < 5  #None 

    while criteria :

        # compute the updated model
        kernel = get_derivative(model)
        model = model_actualisation(model, kernel, step)

        # update condition 
        it += 1
        criteria = it < 5 #None

    return model