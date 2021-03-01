/*
Copyright Kazuya Takahashi 2016-2018.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

module boolean_matrix;

// boolean square matrix
struct BSM(size_t n) if (n > 0)
{
    bool[n][n] s_mat;
    alias s_mat this;

    auto opBinary(string op)(in auto ref BSM!n matrix) const pure
    {
        static if (op == "*")
        {
            BSM!n r;

            foreach (i, ref row; s_mat)
            {
                foreach (j, e; row)
                {
                    if (e)
                    {
                        r[i][] |= matrix[j][];
                    }
                }
            }

            return r;
        }
        else
        {
            static assert(0, `Operator '` ~ op ~ `' is not implemented.`);
        }
    }

    auto opBinary(string op)(in auto ref BCV!n vector) const pure
    {
        static if (op == "*")
        {
            BCV!n r;

            foreach (i, ref row; s_mat)
            {
                foreach (j, e; row)
                {
                    if (e && vector[j])
                    {
                        r[i] = true;
                        break;
                    }
                }
            }

            return r;
        }
        else
        {
            static assert(0, `Operator '` ~ op ~ `' is not implemented.`);
        }
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
    do
    {
        BCV!n v;
        v[k] = true;
        return v;
    }
}

unittest
{
    immutable a = BSM!3([
            [true, false, true], [true, true, false], [false, false, true]
            ]);
    immutable b = BSM!3([
            [false, true, true], [true, false, false], [false, false, true]
            ]);
    immutable c = BSM!3([
            [false, true, true], [true, true, true], [false, false, true]
            ]);

    assert(a * b == c);
}

unittest
{
    immutable a = BSM!3([
            [true, false, true], [true, true, false], [false, false, true]
            ]);
    immutable b = BCV!3([true, false, false]);
    immutable c = BCV!3([true, true, false]);

    assert(a * b == c);
}

unittest
{
    immutable a = BCV!3([false, true, false]);

    assert(BCV!3.e(1) == a);
}

import std.traits : isInstanceOf;

auto transpose(M)(in M a) pure if (isInstanceOf!(BSM, M))
{
    import std.traits : TemplateArgsOf;

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

unittest
{
    immutable a = BSM!3([
            [false, true, false], [false, true, true], [true, false, false]
            ]);
    immutable b = BSM!3([
            [false, false, true], [true, true, false], [false, true, false]
            ]);

    assert(transpose(a) == b);
}

auto isBijective(M)(in M a) pure if (isInstanceOf!(BSM, M))
{
    foreach (i, ref row; a)
    {
        size_t k;

        foreach (j, e; row)
        {
            if (e)
            {
                if (++k == 2)
                {
                    return false;
                }
                else
                {
                    foreach (ref row2; a[i + 1 .. $])
                    {
                        if (row2[j])
                        {
                            return false;
                        }
                    }
                }
            }
        }

        if (k == 0)
        {
            return false;
        }
    }

    return true;
}

unittest
{
    immutable a = BSM!3([
            [true, false, false], [false, false, true], [false, true, false]
            ]);
    immutable b = BSM!3([
            [false, false, true], [false, false, true], [true, true, false]
            ]);
    immutable c = BSM!3([
            [true, false, true], [false, true, false], [false, false, false]
            ]);
    immutable d = BSM!3([
            [false, false, false], [false, true, false], [false, false, false]
            ]);

    assert(isBijective(a));
    assert(!isBijective(b));
    assert(!isBijective(c));
    assert(!isBijective(d));
}

auto isIrreflexive(M)(in M a) pure if (isInstanceOf!(BSM, M))
{
    foreach (i, ref row; a)
    {
        if (row[i])
        {
            return false;
        }
    }

    return true;
}

unittest
{
    immutable a = BSM!3([
            [false, false, true], [false, false, true], [true, true, false]
            ]);
    immutable b = BSM!3([
            [true, false, true], [false, false, true], [true, false, false]
            ]);

    assert(isIrreflexive(a));
    assert(!isIrreflexive(b));
}

auto isSymmetric(M)(in M a) pure if (isInstanceOf!(BSM, M))
{
    foreach (i, ref row; a[0 .. $ - 1])
    {
        foreach (j, e; row[i + 1 .. $])
        {
            if (a[i + j + 1][i] != e)
            {
                return false;
            }
        }
    }

    return true;
}

unittest
{
    immutable a = BSM!3([
            [true, false, true], [false, false, true], [true, true, false]
            ]);
    immutable b = BSM!3([
            [true, false, true], [false, false, true], [true, false, false]
            ]);

    assert(isSymmetric(a));
    assert(!isSymmetric(b));
}
