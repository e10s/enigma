/*
Copyright electrolysis 2016.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

module boolean_matrix;

// boolean square matrix
struct BSM(size_t n) if (n > 0)
{
    bool[n][n] s_mat;
    alias s_mat this;

    auto opMul()(in auto ref BSM!n matrix) const pure
    {
        import std.range : iota, transversal, zip;
        import std.algorithm.searching : canFind;
        import std.algorithm.mutation : copy;
        import std.algorithm.iteration : map;

        BSM!n r;
        foreach (i, ref row; r)
        {
            n.iota.map!(j => zip(s_mat[i][], matrix[].transversal(j))
                .canFind!"a[0]&&a[1]").copy(row[]);
        }

        return r;
    }

    auto opMul()(in auto ref BCV!n vector) const pure
    {
        import std.range : zip;
        import std.algorithm.searching : canFind;
        import std.algorithm.mutation : copy;
        import std.algorithm.iteration : map;

        BCV!n r;
        s_mat[].map!(row => zip(row[], vector[]).canFind!"a[0]&&a[1]")
            .copy(r[]);

        return r;
    }

    static enum I = {
        BSM!n id;
        foreach (i, ref row; id)
        {
            row[i] = true;
        }
        return id;
    }();
}

// boolean column vector
struct BCV(size_t n) if (n > 0)
{
    bool[n] c_vec;
    alias c_vec this;

    static auto e(size_t k) pure
    in
    {
        assert(k < n);
    }
    body
    {
        BCV!n v;
        v[k] = true;
        return v;
    }
}

unittest
{
    immutable a = BSM!3([[true, false, true], [true, true, false],
        [false, false, true]]);
    immutable b = BSM!3([[false, true, true], [true, false, false],
        [false, false, true]]);
    immutable c = BSM!3([[false, true, true], [true, true, true],
        [false, false, true]]);

    assert(a * b == c);
}

unittest
{
    immutable a = BSM!3([[true, false, true], [true, true, false],
        [false, false, true]]);
    immutable b = BCV!3([true, false, false]);
    immutable c = BCV!3([true, true, false]);

    assert(a * b == c);
}

import std.traits : TemplateArgsOf, TemplateOf;

auto transpose(M)(in M a) pure if (__traits(isSame, TemplateOf!M, BSM))
{
    BSM!(TemplateArgsOf!M) ta;
    foreach (i, row; a)
    {
        foreach (j, e; row)
        {
            ta[j][i] = e;
        }
    }
    return ta;
}

auto isSymmetric(M)(in M a) pure if (__traits(isSame, TemplateOf!M, BSM))
{
    foreach (i, row; a[0 .. $-1])
    {
        foreach (j, e; row[1 .. $])
        {
            if (row[j+1] != e)
            {
                return false;
            }
        }
    }

    return true;
}

auto lowerRotator(size_t n)(size_t shift) @property pure
{
    BSM!n r;
    if (shift % n == 0)
    {
        r = BSM!n.I;
    }
    else
    {
        foreach (i, ref e; r)
        {
            e[(i + n - (shift % n)) % n] = true;
        }
    }
    return r;
}

auto upperRotator(size_t n)(size_t shift) @property pure
{
    BSM!n r;
    if (shift % n == 0)
    {
        r = BSM!n.I;
    }
    else
    {
        foreach (i, ref e; r)
        {
            e[(i + (shift % n)) % n] = true;
        }
    }
    return r;
}

auto permutation(size_t N)(in size_t[] substitution)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.range : iota;

    assert(N.iota.isPermutation(substitution));
}
body
{
    BSM!N p;
    foreach (i, s; substitution)
    {
        p[s][i] = true;
    }
    return p;
}
