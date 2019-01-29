# RankClus
A Data Mining Algorithm in the Heterogeneous Information Network.

## Files
* Rankclus

## Input
We have to input a graph structure.It is read by a folder which includes following. 
- X.txt: nodes list of target type (for clustering).
- Y.txt: nodes list of attribute type.
- WXY.csv: Edge between X and Y.
    - [node id (attribute)],[node id (target)],[weight]
- WYY.csv (if necessary) : Edge beteen Y.
    - [node id (attribute)],[node id (target)],[weight]

## Usage

```console
make
./rankclus.out [File Path] [Cluster Number]
```

## Algorithm
The algorithm was based on the paper:`RankClus:integrating clustering with ranking for heterogeneous information network analysis`.
<img width="432" alt="2018-11-01 14 16 21" src="https://user-images.githubusercontent.com/7589567/47833953-c4dc8800-dde0-11e8-9982-486499f5aea6.png">

Referenced by `Table4` in the [paper](http://zuse9-2.se.cuhk.edu.hk/~hcheng/paper/edbt09_ysun.pdf)
