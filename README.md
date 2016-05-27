#limehmm
limehmm is a library for Hidden Markov Models in linear memory. This library was written as part of a third year Computer Science project (the "mini-mémoire") at the Université Libre de Bruxelles. It currently supports : 

- A linear memory implementation for the Forward and Backward procedures (observation likelihood)
- A linear memory implementation for the Viterbi (decode) algorithm
- Two linear memory training algorithms : Viterbi and Baum-Welch
- Discrete distributions
- Silent states
- Mandatory silent begin state / Optional silent end state 
- Fix / free parameters (choose which parameters of the HMM will be trained)

See last section for some examples.

#Dependencies
A c++11 compiler.
#Using the library in your project
##Clone the repo
Firstly, clone and cd into the root directory of the repo :

```
git clone https://github.com/jhellinckx/limehmm.git
cd limehmm
```

##Compile the library
You will then have to compile the source code (to be able to import all the objects file easily don't forget to delete the generated `hmm_test.o` file which contains a `main` function that may conflict with your own) :

````
cd src/
make && rm hmm_test.o
````

##Include
Suppose we have a `foo.cpp` file in which we want to use the library. Include the library as follow : 

```
#include <hmm.hpp>
```

Then, compile it using c++11 standards. You will also need to specify the library object files (without `hmm_test.o` !) and its include path, where `REPO_ROOT` is the path leading to the cloned repo :

```
g++ --std=c++11 foo.cpp -o foo REPO_ROOT/src/*.o -I REPO_ROOT/src/
```


#Examples
##A simple HMM : casino
### Create the states
Create the distributions. Do it by using an initializer list :

```
DiscreteDistribution fair_dist = {{"H", 0.5}, {"T ", 0.5}};
```

Or by using the overloaded `operator[]` :

```
DiscreteDistribution biased_dist;
biased_dist["H"] = 0.75;
biased_dist["T"] = 0.25;
```

Create the states and give them the previously created distributions. The name of each state has to be unique :

```
State fair("fair", fair_dist);
State biased("biased", biased_dist);
```

###Build the HMM
Create a HMM called "casino" and add the states to it :

```
HiddenMarkovModel casino("casino");
casino.add_state(fair);
casino.add_state(biased);
```

Set transitions between states :

```
casino.add_transition(casino.begin(), fair, 0.5);
casino.add_transition(casino.begin(), biased, 0.5);
casino.add_transition(fair, fair, 0.9);
casino.add_transition(fair, biased, 0.1);
casino.add_transition(biased, biased, 0.9);
casino.add_transition(biased, fair, 0.1);
```

Before using any algorithms on this newly created HMM, you will have to call 

```
casino.brew():
```

in order to initialize the model used by the algorithms. Note that calling `brew()` will be necessary each time you modify the HMM.

###Save and load
If you wish to save your HMM on the disk, call (filename and file extension can be set by passing arguments) :

```
casino.save();
```

To load it again, do as follow :

```
HiddenMarkovModel casino_from_disk;
casino_from_disk.load("casino");
```
Let's now use some algorithms on this casino HMM. Don't forget to call `brew()` before using any of those algorithms ! 

###Forward/Backward algorithm
Let the sequence be :

```
std::vector<std::string> sequence = {"T","H","H","T","T","T","H","H"};
```

The Forward and Backward procedures can be set to execute on the entire sequence or only until a given position [in the sequence] : 

```
std::vector<double> init_fwd = casino.forward(sequence, 1); // Iterate once (first position)
std::vector<double> last_fwd = casino.forward(sequence); // Iterate from first position to last position
std::vector<double> init_bwd = casino.backward(sequence, sequence.size()); // Iterate once (last position)
std::vector<double> last_bwd = casino.backward(sequence); // Iterate from last position to first position
```

They are also used to compute the observation probability :

```
double likelihood = casino.likelihood(sequence);
```

###Decode
The `decode` function returns a pair containing the optimal state path and its likelihood :

```
std::vector<std::string> optimal_path = casino.decode(sequence).first;
```

###Training
Currently, the library provides two linear training algorithms : the Viterbi and Baum-Welch training. Viterbi training is the default, in order to use the Baum-Welch algorithm, use the `set_training` method : 

```
casino.set_training(LinearMemoryBaumWelchTraining(nullptr));
```

To check which training algorithm your HMM is currently using, call : 

```
casino.training_type();
```
which returns a string describing the algorithm in use.

Training your HMM is done by calling the `train` method which takes to 5 parameters. Respectively  :

1. The training sequences (Required). Empty sequences are ignored.
2. A transition pseudocount (Optional). Currently only used in the Viterbi training. 
3. The convergence threshold (Optional). When the training iteration improvement is lesser than this threshold, the training stops.
4. The minimum number of iterations (Optional).
5. The maximum number of iterations (Optional). 

The inner HMM model is modified after each training iteration but the HMM directly accessible values are only set once the training is over. The `train` method returns the global improvement. Example with a transition pseudocount set to 1 :

```
std::vector<std::vector<std::string>> training_sequences = {{"T", "H", "H", "T"}, {"T", "H", "H", "T"}, {"T", "H", "H", "T"}};
double improvement = casino.train(training_sequences, 1.0);
```

More details can be found in `hmm.hpp` and `hmm_test.cpp`.