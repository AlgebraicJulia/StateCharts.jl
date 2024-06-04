module NetworkSchemaInterfaces

export SchemaTheory, Open, obname, obnames, nob, nhom, homame, homnames, attrtypename, attrtypenames, attrname, attrnames

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

nob(s::AbstractSchema) = nparts(s, :Ob)
nhom(s::AbstractSchema) = nparts(s, :Hom)
nat(s::AbstractSchema) = nparts(s, :AttrType)
na(s::AbstractSchema) = nparts(s, :Attr)

obname(s::AbstractSchema, ob) = has_subpart(s, :obname) ? subpart(s, ob, :obname) : (1:nob(s))[ob]
obnames(s::AbstractSchema) = map(ob -> obname(s, ob), 1:nob(s))

homname(s::AbstractSchema, hom) = has_subpart(s, :homname) ? subpart(s, hom, :homname) : (1:nhom(s))[hom]
homnames(s::AbstractSchema) = map(hom -> homname(s, hom), 1:nhom(s))

attrtypename(s::AbstractSchema, at) = has_subpart(s, :attrtypename) ? subpart(s, at, :attrtypename) : (1:nat(s))[at]
attrtypenames(s::AbstractSchema) = map(at -> attrtypename(s, at), 1:nat(s))

attrname(s::AbstractSchema, a) = has_subpart(s, :attrname) ? subpart(s, a, :attrname) : (1:na(s))[a]
attrnames(s::AbstractSchema) = map(a -> attrname(s, a), 1:na(s))

# define the open schema as structured cospan
const OpenSchemaObUntyped, OpenSchemaUntyped = OpenACSetTypes(SchemaUntyped, :Ob)
const OpenSchemaOb, OpenSchema = OpenSchemaObUntyped{Symbol}, OpenSchemaUntyped{Symbol}

Open(s::SchemaTheory) = OpenSchema(s, map(x -> FinFunction([x], nob(s)), 1:nob(s))...)
Open(s::SchemaTheory, legs...) = begin
  s_idx = Dict(obname(s, ob) => ob for ob in 1:nob(s))
  OpenSchema(s, map(l -> FinFunction(map(i -> s_idx[i], l), nob(s)), legs)...)
end
Open(n, s::SchemaTheory, m) = Open(s,n,m)

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

end