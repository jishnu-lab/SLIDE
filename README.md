## Contents

- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [Hardware requirements](#Hardware requirements)







## Installation Procedure
# SLIDE
First install `SLIDE` in Python by running the following code from a  R command terminal:



```library(devtools)```   
```install_github("jishnu-lab/SLIDE")```


To run the slide function, use the following command:

```library(doParallel)```

For capturing the interaction effect ```do_interacts = TRUE```

## Demo
From the  ```./test``` folder read the ```z``` and ```y``` using this command

```z <- ReadRDS('./test/z.rds')```

```y <- ReadRDS('./test/y.rds')``` 

```res <- SLIDE(z,y,method = 4,do_interacts = TRUE)```

## Results
The expected outcome for the run should be:

**For marginal variables:**

```res$marginal_vars```

```"z4"  "z10" "z46" ```

**For the interactors:**

```Z4.Z3     Z4.Z12     Z4.Z18       Z4.Z27     Z10.Z18     Z46.Z41 ```



## Hardware requirements
For minimal performance of the ```SLIDE```,  a computer with about 8 GB of RAM is required. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The expected run time with computer with 32GB ram and 16 cores for the demo command to run is 

```282.25 seconds ```




   
