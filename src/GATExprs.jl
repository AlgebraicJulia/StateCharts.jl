module GATExprs

import Base: +, -, *, exp
import Catlab: ⋅
using Catlab

@theory ThStateChartsAlgebra begin
    State::TYPE
    Quantity::TYPE

    # Multiplicative commutative monoid
    @op (⋅) := times
    times(x::Quantity, y::Quantity)::Quantity
    unit()::Quantity

    x ⋅ (y ⋅ z) == (x ⋅ y) ⋅ z ⊣ [(x,y,z)::Quantity]
    x ⋅ unit() == x ⊣ [x::Quantity]
    unit() ⋅ x == x ⊣ [x::Quantity]
    x ⋅ y == y ⋅ x ⊣ [(x,y)::Quantity]

    # Additive abelian group
    @op (+) := plus
    plus(x::Quantity, y::Quantity)::Quantity
    zero()::Quantity
    @op (~) := neg
    neg(x::Quantity)::Quantity

    x + (y + z) == (x + y) + z ⊣ [(x,y,z)::Quantity]    
    x + zero() == x ⊣ [x::Quantity]
    zero() + x == x ⊣ [x::Quantity]
    x + y == y + x ⊣ [(x,y)::Quantity] 
    x + neg(x) == zero() ⊣ [(x,y)::Quantity]
    neg(x) + x == zero() ⊣ [(x,y)::Quantity]
    # distributive law 
    x ⋅ (y + z) == x ⋅ y + x ⋅ z ⊣ [(x,y,z)::Quantity]

    # Convenience term for subtraction
    @op (-) := minus
    minus(x::Quantity, y::Quantity)::Quantity
    x - y == x + (~y) ⊣ [(x,y)::Quantity]
    x - x == zero() ⊣ [x::Quantity]

    # exponential 
    exp(x::Quantity)::Quantity




end





end