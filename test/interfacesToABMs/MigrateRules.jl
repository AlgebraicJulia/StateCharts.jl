module TestMigrateRules

using Test
using StateCharts
using Catlab

# Create a schema with S, I and R as objects
@present SchSIR(FreeSchema) begin
    Person::Ob
    (S,I,R)::Ob
    SPerson::Hom(S,Person)
    IPerson::Hom(I,Person)
    RPerson::Hom(R,Person)
end

SchSIRACSet=schemaACSet(SchSIR)

@test  obnames(SchSIRACSet) == [:Person, :S,:I,:R]
@test  homnames(SchSIRACSet) == [:SPerson,:IPerson,:RPerson]

# Create a schema with a state chart as an attribute
@present SchSIR2(FreeSchema) begin
    Person::Ob
    
    Inf::AttrType

    inf::Attr(Person,Inf)
end

SchSIRACSet2=schemaACSet(SchSIR2)

@test obnames(SchSIRACSet2) == [:Person]
@test attrtypenames(SchSIRACSet2) == [:Inf]
@test attrnames(SchSIRACSet2) == [:inf]


# test basic schema convert to acset
basicschema = BasicSchema([:P,:V],[(:PV,:P,:V)],[:Inf,:Coord],[(:inf,:P,:Inf),(:coord,:V,:Coord)])
SchACSetBS = schemaACSet(basicschema)

@test obnames(SchACSetBS) == [:P,:V]
@test homnames(SchACSetBS) == [:PV]
@test attrtypenames(SchACSetBS) == [:Inf,:Coord]
@test attrnames(SchACSetBS) == [:inf,:coord]

# test convert back to present schema
@test SchSIR == schemaPresent(schemaACSet(SchSIR))
@test SchSIR2 == schemaPresent(schemaACSet(SchSIR2))
end