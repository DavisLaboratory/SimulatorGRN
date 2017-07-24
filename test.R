source('simulatorS4-methods.R')
source('simulatorS4.R')

#----case 1----
n1=new('Node', name='A')
n2=new('Node', name='B')
n3=new('Node', name='C')
n4=new('Node', name='D')
e1=new('EdgeOr', from=c(n2), to=n1, activation=F)
e2=new('EdgeOr', from=c(n3), to=n1)
e3=new('EdgeOr', from=c(n4), to=n1)
e4=new('EdgeOr', from=c(n3), to=n2)
and1=new('EdgeAnd', from=c(n1,n2,n3), to=n4)

n1$inedges=list(e1,e2,e3)
n2$inedges=list(e4)
n4$inedges=list(and1)

print(generateRateEqn(n1))
print(generateRateEqn(n2))
print(generateRateEqn(n3))
print(generateRateEqn(n4))

#----case 2----
grn=new('GraphGRN')
grn=addNode(grn, 'A')
grn=addNode(grn, 'B')
grn=addNode(grn, 'C')
grn=addNode(grn, 'D')
grn=addEdge(grn,'B','A',activation=F)
grn=addEdge(grn,'C','A')
grn=addEdge(grn,'D','A')
grn=addEdge(grn,'C','B')
grn=addEdge(grn,c('A','B','C'),'D',edgetype='and')

print(generateRateEqn(grn@nodeset[['A']]))
print(generateRateEqn(grn@nodeset[['B']]))
print(generateRateEqn(grn@nodeset[['C']]))
print(generateRateEqn(grn@nodeset[['D']]))

