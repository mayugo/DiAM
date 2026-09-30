import numpy as np
sigma = np.array([[900, 200, 0], 
                  [200, 300, 0],
                  [  0,   0, 0]])   # Estat de tensions
p, v = np.linalg.eig(sigma)         # Calcular valors propis 'p' i vectors propis 'v'