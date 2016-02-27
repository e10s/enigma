/*
Copyright electrolysis 2016.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

module meta_workaround;

static if (is(typeof({ import std.meta : Repeat;  })))
{
    static import std.meta;

    ///
    alias Repeat = std.meta.Repeat;
}
else
{
    /// From Phobos pre-release.
    template Repeat(size_t n, TList...) if (n > 0)
    {
        import std.meta : AliasSeq;

        static if (n == 1)
        {
            alias Repeat = AliasSeq!TList;
        }
        else static if (n == 2)
        {
            alias Repeat = AliasSeq!(TList, TList);
        }
        else
        {
            alias R = Repeat!((n - 1) / 2, TList);
            static if ((n - 1) % 2 == 0)
            {
                alias Repeat = AliasSeq!(TList, R, R);
            }
            else
            {
                alias Repeat = AliasSeq!(TList, TList, R, R);
            }
        }
    }
}
