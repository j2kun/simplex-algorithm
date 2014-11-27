from test import test

import simplex

def testFromPost():
   cost = [1,1,1]
   gts = [[0,1,4]]
   gtB = [10]
   lts = [[3,-2,0]]
   ltB = [7]
   eqs = [[1,1,0]]
   eqB = [2]

   expectedCost = [1,1,1,0,0]
   expectedConstraints = [[0,1,4,-1,0], [3,-2,0,0,1], [1,1,0,0,0]]
   expectedThresholds = [10,7,2]
   test((expectedCost, expectedConstraints, expectedThresholds),
         simplex.standardForm(cost, gts, gtB, lts, ltB, eqs, eqB))


def test2():
   cost = [1,1,1]
   lts = [[3,-2,0]]
   ltB = [7]
   eqs = [[1,1,0]]
   eqB = [2]

   expectedCost = [1,1,1,0]
   expectedConstraints = [[3,-2,0,1], [1,1,0,0]]
   expectedThresholds = [7,2]
   test((expectedCost, expectedConstraints, expectedThresholds),
         simplex.standardForm(cost, lessThans=lts, ltThreshold=ltB, equalities=eqs, eqThreshold=eqB))


def test3():
   cost = [1,1,1]
   eqs = [[1,1,0], [2,2,2]]
   eqB = [2, 5]

   expectedCost = [1,1,1]
   expectedConstraints = [[3,-2,0], [1,1,0]]
   expectedThresholds = [2, 5]
   test((expectedCost, expectedConstraints, expectedThresholds),
         simplex.standardForm(cost, equalities=eqs, eqThreshold=eqB))


if __name__ == "__main__":
   testFromPost()
