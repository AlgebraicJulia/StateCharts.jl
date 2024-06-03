module NetworkSchemaInterfaces

using Test
using StateCharts
using Catlab

# test schema of schema
ob=(:S,:I,:R,:V)
hom=(:SV=>(:S,:V),
     :IV=>(:I,:V),
     :RV=>(:R,:V))


schemaCSet=SchemaTheory(ob, hom, [], [])

ob2=(:V)
hom2=[]
attrtype=[:SIR]
attr=[:VSIR=>(:V,:SIR)]

schemaACSet=SchemaTheory(ob2, hom2, attrtype,attr)

# test composition of schemas - example 1: object person (which is V in this example) in state chart generated schema is identified
# with the vertice object in basic graph.

# Define basic graph schema
BasicGraphSchema=SchemaTheory([:E,:V], [:src=>(:E,:V),:tgt=>(:E,:V)], [], [])

# Compose the two schemas
composedSchema=apex(compose(Open([:S],schemaCSet,[:V]),Open([:V],BasicGraphSchema,[:E])))

@test  obnames(composedSchema) == [:S,:I,:R,:V,:E]
@test  homnames(composedSchema) == [:SV,:IV,:RV,:src,:tgt]

# test composition of schemas - example 2: there is a morphism from the object person in state chart generated schema to
# the vertice object in basic graph.

schemaSIRP = SchemaTheory([:S,:I,:R,:P], [:SP=>(:S,:P),:IP=>(:I,:P),:RP=>(:R,:P)], [], [])
schemaPV = SchemaTheory([:P,:V], [:PV=>(:P,:V)], [], [])
composedSchema2 = (Open([:S],schemaSIRP,[:P]) · Open([:P],schemaPV,[:V])) · Open([:V],BasicGraphSchema,[:E]) |> apex
# the syntax of compose can also using julia unicode of \cdotp

@test  obnames(composedSchema2) == [:S,:I,:R,:P,:V,:E]
@test  homnames(composedSchema2) == [:SP,:IP,:RP,:PV,:src,:tgt]

end