# Sequence Space

```{warning}
This library is under developement
```

This module solves dynamic systems using the sequence space approach. 

The goal is to find: 
$$F(X_t,\mathcal{E}_t)=0$$

Usually $F$ is composed by:

$$f_t(x_{t-1},x_{t},x_{t+1},\epsilon_t)$$

This simplifies the problem since the dependence of $x_t$ to $x_{t+2}$ is only though $x_{t+1}$. By creating new variables  we can acomodate further leads and lags. 

TODO: the function only depends on $\epsilon_t$ or also $\epsilon_{t-1}$ and $\epsilon_{t+1}$?

```{eval-rst}
.. autofunction:: sequence_space.sequence_space.build_F
```


```{eval-rst}
.. autofunction:: sequence_space.sequence_space.solve_model
```

