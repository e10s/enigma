/*
Copyright Kazuya Takahashi 2018.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

module structure;

struct SetElement(size_t cardinalityOfSet) if (cardinalityOfSet > 0)
{
    import boolean_matrix : BCV;

    immutable size_t id;

    this(size_t id)
    {
        this.id = id;
    }

    this(BCV!cardinalityOfSet entity)
    {
        import std.algorithm.searching : countUntil;

        id = entity[].countUntil(true);
    }

    private @property entity() const
    {
        return cast(immutable) BCV!cardinalityOfSet.e(id);
    }
}

unittest
{
    import boolean_matrix : BCV;

    auto elem1 = SetElement!4(2);
    assert(elem1.id == 2);
    assert(elem1.entity == BCV!4([false, false, true, false]));

    auto elem2 = SetElement!5(BCV!5([false, false, false, true, false]));
    assert(elem2.id == 3);
    assert(elem2.entity == BCV!5([false, false, false, true, false]));
}

struct PermutationElement(size_t degree) if (degree > 0)
{
    import boolean_matrix : BSM;

    immutable(BSM!degree) entity;
    this(BSM!degree entity)
    {
        this.entity = entity;
    }

    auto opBinary(string op)(in auto ref PermutationElement!degree g) const pure
    {
        static if (op == "*")
        {
            return immutable(PermutationElement!degree)(this.entity * g.entity);
        }
        else
        {
            static assert(0, `Operator '` ~ op ~ `' is not implemented.`);
        }
    }

    auto opBinary(string op)(in auto ref SetElement!degree s) const pure
    {
        static if (op == "*")
        {
            return immutable(SetElement!degree)(this.entity * s.entity);
        }
        else
        {
            static assert(0, `Operator '` ~ op ~ `' is not implemented.`);
        }
    }

    auto transpose() const pure
    {
        import boolean_matrix : transpose;

        return immutable(PermutationElement!degree)(entity.transpose);
    }

    // Delegates "isXXX" properties to boolean_matrix.
    @property opDispatch(string s)() const pure if (s.length > 2 && s[0 .. 2] == "is")
    {
        import boolean_matrix;

        return mixin(s)(entity);
    }
}

unittest
{
    import boolean_matrix : BSM;

    immutable a = PermutationElement!3(BSM!3([
                [true, false, true], [true, true, false], [false, false, true]
            ]));
    immutable b = PermutationElement!3(BSM!3([
                [false, true, true], [true, false, false], [false, false, true]
            ]));
    immutable c = PermutationElement!3(BSM!3([
                [false, true, true], [true, true, true], [false, false, true]
            ]));

    assert(a * b == c);
}

unittest
{
    import boolean_matrix : BCV, BSM;

    immutable a = PermutationElement!3(BSM!3([
                [true, false, true], [true, true, false], [false, false, true]
            ]));
    immutable b = SetElement!3(BCV!3([true, false, false]));
    immutable c = SetElement!3(BCV!3([true, true, false]));

    assert(a * b == c);
}

auto cyclicPermutation(size_t degree)(size_t shift) pure
{
    import boolean_matrix : BSM;

    BSM!degree r;
    if (shift % degree == 0)
    {
        r = BSM!degree.I;
    }
    else
    {
        foreach (i, ref e; r)
        {
            e[(i + (shift % degree)) % degree] = true;
        }
    }
    return PermutationElement!degree(r);
}

unittest
{
    import boolean_matrix : BSM;

    immutable a = PermutationElement!3(BSM!3([
                [true, false, false], [false, true, false], [false, false, true]
            ]));
    immutable b = PermutationElement!3(BSM!3([
                [false, true, false], [false, false, true], [true, false, false]
            ]));

    assert(cyclicPermutation!3(0) == a);
    assert(cyclicPermutation!3(21) == a);
    assert(cyclicPermutation!3(1) == b);
    assert(cyclicPermutation!3(22) == b);
}

auto cyclicPermutationInv(size_t degree)(size_t shift) pure
{
    import boolean_matrix : BSM;

    BSM!degree r;
    if (shift % degree == 0)
    {
        r = BSM!degree.I;
    }
    else
    {
        foreach (i, ref e; r)
        {
            e[(i + degree - (shift % degree)) % degree] = true;
        }
    }
    return PermutationElement!degree(r);
}

unittest
{
    import boolean_matrix : BSM;

    immutable a = PermutationElement!3(BSM!3([
                [true, false, false], [false, true, false], [false, false, true]
            ]));
    immutable b = PermutationElement!3(BSM!3([
                [false, false, true], [true, false, false], [false, true, false]
            ]));

    assert(cyclicPermutationInv!3(0) == a);
    assert(cyclicPermutationInv!3(21) == a);
    assert(cyclicPermutationInv!3(1) == b);
    assert(cyclicPermutationInv!3(22) == b);
}

unittest
{
    assert(cyclicPermutation!6(3) == cyclicPermutationInv!6(3).transpose);
    assert(cyclicPermutation!6(51) == cyclicPermutationInv!6(51).transpose);
    assert(cyclicPermutation!6(1).transpose == cyclicPermutationInv!6(1));
    assert(cyclicPermutation!6(22).transpose == cyclicPermutationInv!6(22));
}

auto permutation(size_t degree)(in size_t[] substitution)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.range : iota;

    assert(degree.iota.isPermutation(substitution));
}
do
{
    import boolean_matrix : BSM;

    BSM!degree p;
    foreach (i, s; substitution)
    {
        p[s][i] = true;
    }
    return PermutationElement!degree(p);
}

unittest
{
    import boolean_matrix : BSM;

    immutable a = PermutationElement!3(BSM!3([
                [true, false, false], [false, true, false], [false, false, true]
            ]));
    immutable b = PermutationElement!3(BSM!3([
                [false, true, false], [false, false, true], [true, false, false]
            ]));
    immutable c = PermutationElement!3(BSM!3([
                [false, false, true], [false, true, false], [true, false, false]
            ]));

    assert(permutation!3([0, 1, 2]) == a);
    assert(permutation!3([2, 0, 1]) == b);
    assert(permutation!3([2, 1, 0]) == c);
}
