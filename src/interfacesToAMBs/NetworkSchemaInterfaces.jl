module NetworkSchemaInterfaces

using Catlab

vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n::SubArray) = collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

######################## define the schema of schema to combine multiple schema ################
@present TheorySchema(FreeSchema) begin
  Ob::Ob
  Hom::Ob
  Attr::Ob
  AttrType::Ob

  homsrc::Hom(Hom,Ob)
  homtgt::Hom(Hom,Ob)
  attrsrc::Hom(Attr,Ob)
  attrtgt::Hom(Attr,AttrType)

# Attributes:
  Name::AttrType
  
  obname::Attr(Ob, Name)
  homname::Attr(Hom, Name)
  attrname::Attr(Attr, Name)
  attrtypename::Attr(AttrType, Name)
end

@abstract_acset_type AbstractSchema
@acset_type SchemaUntyped(TheorySchema, index=[:homsrc,:homtgt,:attrsrc,:attrtgt]) <: AbstractSchema
const SchemaTheory = SchemaUntyped{Symbol} 

function SchemaTheory(ob,hom,attrtype,attr)
    st = SchemaTheory()

    ob = vectorify(ob)
    hom = vectorify(hom)
    attrtype = vectorify(attrtype)
    attr = vectorify(attr)

    add_parts!(st, :Ob, length(ob); obname = ob)
    add_parts!(st, :AttrType, length(attrtype); attrtypename = attrtype)

    ob_idx=state_dict(ob)
    atype_idx=state_dict(attrtype)

    hom_name=map(first,hom)
    hom_ob=map(last,hom)
    hsrc = map(x->ob_idx[x], map(first, hom_ob))
    htgt = map(x->ob_idx[x], map(last, hom_ob))
    add_parts!(st, :Hom, length(hom); homsrc=hsrc, homtgt=htgt,homname=hom_name)

    attr_name=map(first,attr)
    attr_oba=map(last,attr)
    asrc = map(x->ob_idx[x], map(first, attr_oba))
    atgt = map(x->atype_idx[x], map(last, attr_oba))
    add_parts!(st, :Attr, length(attr); attrsrc=asrc, attrtgt=atgt,attrname=attr_name)

    return st
end



####### example of 
ob=(:S,:I,:R,:V)
hom=(:SV=>(:S,:V),
     :IV=>(:I,:V),
     :RV=>(:R,:V))

ob2=(:V)
hom2=[]
attrtype=[:SIR]
attr=[:VSIR=>(:V,:SIR)]

a=SchemaTheory(ob, hom, [], [])

b=SchemaTheory(ob2, hom2, attrtype,attr)


  
end