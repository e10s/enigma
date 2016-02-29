/*
Copyright Kazuya Takahashi 2016.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

void main()
{
    import std.stdio;
    import enigma;

    auto m3 = EnigmaM3(entryWheelABC, rotorIII, rotorII, rotorI, reflectorB, "AAA");
    foreach (c; readln!dstring)
    {
        write(m3(c));
    }
}
