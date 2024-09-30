# this module defines some commonly used migration rules
module MigrateRules

export schemaACSet, schemaPresent

using Catlab
using ..NetworkSchemaInterfaces

# schema convert functions between ACSet and schemas
# create an acset schema from a present schema
function schemaACSet(s::Presentation)
    obs = nameof.(generators(s,:Ob))
    h = generators(s,:Hom)
    homs = map((name,dom,codom)->name=>(dom,codom),nameof.(h),nameof.(dom.(h)),nameof.(codom.(h)))
    attrtypes = nameof.(generators(s,:AttrType))
    a = generators(s,:Attr)
    attrs = map((name,dom,codom)->name=>(dom,codom),nameof.(a),nameof.(dom.(a)),nameof.(codom.(a)))
    schema = SchemaTheory(obs,homs,attrtypes,attrs)
    schema
end

# create an acset schema from a basic schema
function schemaACSet(s::BasicSchema)
    obs = objects(s)
    hs = [a=>(b,c) for (a,b,c) in homs(s)]
    ats = attrtypes(s)
    as = [a=>(b,c) for (a,b,c) in attrs(s)]
    schema = SchemaTheory(obs,hs,ats,as)
    schema
end

# create a @present schema from ACSet schema
function schemaPresent(s)
    pres = Presentation(FreeSchema)
    for ob in parts(s, :Ob)
        add_generator!(pres, Ob(FreeSchema.Ob, s[ob, :obname]))
    end
    obs = generators(pres, :Ob)
    for hom in parts(s, :Hom)
        f = Hom(s[hom, :homname], obs[s[hom, :homsrc]], obs[s[hom, :homtgt]])
        add_generator!(pres, f)
    end
    for attrtype in parts(s, :AttrType)
        add_generator!(pres, AttrType(FreeSchema.AttrType, s[attrtype, :attrtypename]))
    end
    attrtypes = generators(pres, :AttrType)
    for attr in parts(s, :Attr)
        f = Attr(s[attr, :attrname], obs[s[attr, :attrsrc]], attrtypes[s[attr, :attrtgt]])
        add_generator!(pres, f)
    end
    pres
end

end