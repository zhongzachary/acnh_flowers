# ACNH Flowers Cross Breeding Calculator

This program allows you to plan how you want to cross breed flowers in Animal Crossing New Horizon.
It provides the probability for your cross breeding results.

The genetic information comes from Paleh ([reddit](https://www.reddit.com/r/AnimalCrossing/comments/fsyn1x/acnhacnl_flower_genetics_guide/), [Google Docs](https://docs.google.com/document/d/1ARIQCUc5YVEd01D7jtJT9EEJF45m07NXhAm4fOpNvCs)).

## Example Usages
I recommend you to use repl.it to run it at https://acnhflowers.zhongzachary.repl.run.

Once you fire up the online python console via the link above (and wait for `pandas` to get installed), you can run the following. 

(Alternatively, you can download the package at my [Github](https://github.com/zhongzachary/acnh_flowers) and run `python -i acnh_flowers.py`.)

The following flower types are supported.
```
Cosmo
Hyacinth
Lily
Mum
Pansy
Rose
Tulip
Windflower
```
#### See Gene Sequences of A Flower Type
```
>>> Table(Mum)
```
Output
```
Red   |Aa Aa Aa
... (multiple rows omitted)
Purple|aa aa AA
      |Aa AA aa
      |Aa AA Aa
      |Aa AA AA
      |AA Aa aa
      |AA Aa Aa
Green |AA AA aa
      |AA AA Aa
```

#### Start Breeding
```
>>> plan = Breeding(Mum)
```
This initiates the breeding plan for Mum. 
Note that each id represents a breeding node (which I will explain later). 
All 3 seeds are included.
- `gen` (generation) is 0 for seeds
- they have no `parents`
- `prob` (probability of getting this) is 1 (since you can always buy it)
- `gsp` (gene sequence probability) is the unique gene for mum's seeds with probability 1 for each sequence.
```
id   |color |species |prob   |rarity |gen|parents|gsp
    0|Red   |Mum     |1.00000|      0|  0|       |AA aa aa   1.0000
    1|Yellow|Mum     |1.00000|      0|  0|       |aa AA aa   1.0000
    2|White |Mum     |1.00000|      0|  0|       |aa aa Aa   1.0000
```
#### Choose What to Breed Next
Say you want to breed the red and yellow seeds.
```
>>> plan.breed(0,1)  # breeding node 0 and node 1
```
Output
```
id   |color |species |prob   |rarity |gen|parents|gsp
    0|Red   |Mum     |1.00000|      0|  0|       |AA aa aa   1.0000
    1|Yellow|Mum     |1.00000|      0|  0|       |aa AA aa   1.0000
    2|White |Mum     |1.00000|      0|  0|       |aa aa Aa   1.0000
    3|Yellow|Mum     |1.00000|    0.0|  1|    1,0|Aa Aa aa   1.0000
```
Note that 3 is created. 
This tells you yellow flower with that gene is the guaranteed output.

Let's breed 3 with itself.
```
>>> plan.breed(3,3)
```
Output 
```
id   |color |species |prob   |rarity |gen|parents|gsp
    0|Red   |Mum     |1.00000|      0|  0|       |AA aa aa   1.0000
    1|Yellow|Mum     |1.00000|      0|  0|       |aa AA aa   1.0000
    2|White |Mum     |1.00000|      0|  0|       |aa aa Aa   1.0000
    3|Yellow|Mum     |1.00000|    0.0|  1|    1,0|Aa Aa aa   1.0000
    4|White |Mum     |0.06250|2.77258|  2|      3|aa aa aa   1.0000
    5|Pink  |Mum     |0.12500|2.07944|  2|      3|Aa aa aa   1.0000
    6|Red   |Mum     |0.06250|2.77258|  2|      3|AA aa aa   1.0000
    7|Yellow|Mum     |0.43750|0.82667|  2|      3|aa Aa aa   0.2857
                                                  Aa Aa aa   0.5714
                                                  aa AA aa   0.1429
    8|Purple|Mum     |0.25000|1.38629|  2|      3|AA Aa aa   0.5000
                                                  Aa AA aa   0.5000
    9|Green |Mum     |0.06250|2.77258|  2|      3|AA AA aa   1.0000
```
This tells us the gen-1 yellow mum, if cross bred with itself, can create 6 different colors of mums (4 to 9).
It has a 0.0625 (6.25%) probability of getting a green mum; 0.25 (25%) probability of getting a purple mum.
Note that the `gsp` column tells you the gene composition of that breeding node.
For example, half of gen-2 purple flowers (breeding node 8) are `AA Aa aa` and the other half are `Aa AA aa`.

## Technical Notes
#### Breeding Nodes
Each id in the output represents a breeding node.
It represents all the flowers from the same parents and having the same colors.
Flowers within the same node can have different gene sequences.

### `breed` Function
When you are using the `breed(id1, id2)` function to breed two nodes, it adds all the possible children nodes to the plan.
However, the following will not be added:
- any children nodes that have the same gene sequence probability table as either of the parents
- if you breed the same parents twice

#### Gene Sequences
Instead of representing gene sequences using different letters (Rr, Yy, etc.), I use Aa here and use the location of genes to differentiate the type of genes.

