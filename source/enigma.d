// Written in the D programming language.

/*
Copyright Kazuya Takahashi 2016.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*/
/**
 * A library for simulating the Enigma machines.
 *
 * Copyright: Copyright Kazuya Takahashi 2016.
 * License:   $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0).
 * Authors:   Kazuya Takahashi
 */
module enigma;

private enum size_t N = 26;

private template isSomeStringOrDcharRange(T)
{
    import std.range.primitives : isInfinite, isInputRange, ElementType;
    import std.traits : isSomeString;

    enum isSomeStringOrDcharRange = isSomeString!T ||
        (isInputRange!T && !isInfinite!T && is(ElementType!T : dchar));
}


///
struct Rotor
{
    import boolean_matrix : BSM;

    immutable BSM!N perm;
    private immutable bool hasNotch = false;
    private immutable size_t[] turnovers;

    /++
     + Constructs a rotor with no turnover notches.
     + If ringOffset is `2`, it corresponds to "C-03".
     +/
    this()(in auto ref BSM!N perm, size_t ringOffset) pure
    {
        import boolean_matrix : lowerRotator, upperRotator;

        this.perm = lowerRotator!N(ringOffset) * perm * upperRotator!N(ringOffset);
    }

    /++
     + Constructs a rotor with one turnover notch.
     + If turnover is `1`, the next rotor steps when this rotor steps from B to C.
     + If ringOffset is `2`, it corresponds to "C-03".
     +/
    this()(in auto ref BSM!N perm, size_t turnover, size_t ringOffset) pure
    {
        this(perm, ringOffset);
        this.turnovers = [turnover % N];
        hasNotch = true;
    }

    /++
     + Constructs a rotor with two turnover notches.
     + If turnover1 is `1` and turnover2 is `25`, the next rotor steps
     + when this rotor steps from B to C and from Z to A.
     + If ringOffset is `2`, it corresponds to "C-03".
     +/
    this()(in auto ref BSM!N perm, size_t turnover1, size_t turnover2, size_t ringOffset) pure
    {
        import std.algorithm.sorting : sort;
        import std.array : array;

        this(perm, ringOffset);
        this.turnovers = [turnover1 % N, turnover2 % N].sort().array;
        hasNotch = true;
    }
}

/++
 + A convenience function to make a rotor with no turnover notches
 + from a forward substitution.
 + If ringSetting is `'C'`, it corresponds to "C-03".
 +/
auto rotor(S)(in S forwardSubstitution, dchar ringSetting) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength;

    assert(forwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(forwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    assert(ringSetting.isAlpha);
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Rotor(forwardSubstitution.map!toUpper.map!"a-'A'".array.permutation!N,
        ringSetting.toUpper - 'A');
}

/++
 + A convenience function to make a rotor with one turnover notch
 + from a forward substitution.
 + If turnover is `'B'`, the next rotor steps when this rotor steps from B to C.
 + If ringSetting is `'C'`, it corresponds to "C-03".
 +/
auto rotor(S)(in S forwardSubstitution, dchar turnover, dchar ringSetting) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength;

    assert(forwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(forwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    assert(turnover.isAlpha);
    assert(ringSetting.isAlpha);
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Rotor(forwardSubstitution.map!toUpper.map!"a-'A'".array.permutation!N,
        turnover.toUpper - 'A', ringSetting.toUpper - 'A');
}

/++
 + A convenience function to make a rotor with two turnover notches
 + from a forward substitution.
 + If turnover1 is `'B'` and turnover2 is `'Z'`, the next rotor steps
 + when this rotor steps from B to C and from Z to A.
+ If ringSetting is `'C'`, it corresponds to "C-03".
 +/
auto rotor(S)(in S forwardSubstitution, dchar turnover1, dchar turnover2, dchar ringSetting) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength;

    assert(forwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(forwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    assert(turnover1.isAlpha);
    assert(turnover2.isAlpha);
    assert(ringSetting.isAlpha);
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Rotor(forwardSubstitution.map!toUpper.map!"a-'A'".array.permutation!N,
        turnover1.toUpper - 'A', turnover2.toUpper - 'A', ringSetting.toUpper - 'A');
}

///
struct Plugboard
{
    import boolean_matrix : BSM;

    immutable BSM!N perm;
    ///
    this()(in auto ref BSM!N perm) pure
    {
        this.perm = cast(immutable) perm;
    }

    alias perm this;
}

/++
 + A convenience function to make a plugboard from a substitution.
 +/
auto plugboard(S)(in S substitution) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : toUpper;
    import std.range : iota, walkLength;

    assert(substitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(substitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
}
out (r)
{
    import boolean_matrix : transpose;

    assert(r == r.transpose);
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Plugboard(substitution.map!toUpper.map!"a-'A'".array.permutation!N);
}

///
struct Reflector
{
    import boolean_matrix : BSM;

    immutable BSM!N perm;
    ///
    this()(in auto ref BSM!N perm) pure
    {
        this.perm = cast(immutable) perm;
    }

    alias perm this;
}

/++
 + A convenience function to make a reflector from a substitution.
 +/
auto reflector(S)(in S substitution) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : toUpper;
    import std.range : iota, walkLength;

    assert(substitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(substitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
}
out (r)
{
    import boolean_matrix : transpose;

    assert(r == r.transpose);
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Reflector(substitution.map!toUpper.map!"a-'A'".array.permutation!N);
}

/// Currently only M3- or M4-like machines are available.
struct Enigma(size_t rotorN, bool fixedFinalRotor = false)
{
    private immutable Plugboard plugboard;
    private immutable Rotor[rotorN] rotors;
    private immutable Reflector reflector;
    private size_t[rotorN] rotationStates;

    import meta_workaround : Repeat;

    ///
    this(in Plugboard plugboard, in Repeat!(rotorN, Rotor) rotors,
        in Reflector reflector, in dchar[rotorN] rotorStartPos)
    in
    {
        foreach (dchar c; rotorStartPos)
        {
            import std.ascii : isAlpha;

            assert(c.isAlpha);
        }
    }
    body
    {
        foreach (i, ref e; rotationStates)
        {
            import std.ascii : toUpper;

            e = rotorStartPos[i].toUpper - 'A';
        }

        this.plugboard = cast(immutable) plugboard;
        this.rotors[] = cast(immutable)[rotors][];
        this.reflector = cast(immutable) reflector;
    }

    private void step()
    {
        enum movableRotorN = fixedFinalRotor ? rotorN - 1 : rotorN;
        bool[movableRotorN] stepFlag;

        stepFlag[0] = true;

        // Handles double stepping
        foreach (rotorID; 0 .. movableRotorN - 1)
        {
            import std.algorithm.searching : canFind;

            if (rotors[rotorID].turnovers.canFind(rotationStates[rotorID]))
            {
                stepFlag[rotorID] = true;
                stepFlag[rotorID + 1] = true;
            }
        }

        foreach (rotorID, e; stepFlag)
        {
            if (e)
            {
                rotationStates[rotorID] = (rotationStates[rotorID] + 1) % N;
            }
        }
    }

    import boolean_matrix : BSM;

    /+
     + fwdPerm = (Um*Rm*Lm)*(Um-1*Rm-1*Lm-1)*...*(U0*R0*L0)*P
     +         = Um*(Rm*Lm*Um-1)*(Rm-1*Lm-1*Um-2)*...*(R0*L0)*P
     +         = Um*(Rm*RELm)*(Rm-1*RELm-1)*...*(R0*L0)*P
     +/
    private BSM!N composeForwardPermutation(in ref BSM!N prevPerm, size_t rotorID)
    {
        import boolean_matrix : lowerRotator, upperRotator, Identity;

        ptrdiff_t x = rotorID == 0 ? rotationStates[0] : rotationStates[rotorID] - rotationStates[rotorID - 1];
        BSM!N relRotator = x > 0 ? lowerRotator!N(x) : upperRotator!N(-x);

        auto composedPerm = rotors[rotorID].perm * relRotator * prevPerm;
        return rotorID == rotorN - 1 ? upperRotator!N(rotationStates[rotorID]) * composedPerm
            : composeForwardPermutation(composedPerm, rotorID + 1);
    }

    private auto process(size_t keyInputID)
    {
        step();

        import boolean_matrix : transpose, BCV;

        auto fwdPerm = composeForwardPermutation(plugboard, 0);
        // bwdPerm = fwdPerm^-1 = fwdPerm^T
        auto wholePerm = fwdPerm.transpose * reflector * fwdPerm;
        BCV!N v;
        v[keyInputID] = true;
        auto w = wholePerm * v;
        import std.algorithm.searching : countUntil;

        auto r = w[].countUntil(true);
        assert(r >= 0);
        return r;

    }

    /// Enciphers only an alphabetical character through the current Enigma machine.
    dchar opCall(dchar keyInput)
    {
        import std.ascii : isAlpha, toUpper;

        return keyInput.isAlpha ? process(keyInput.toUpper - 'A') + 'A' : keyInput;
    }
}

/// The M3 Enigma machine, which has three rotor slots.
alias EnigmaM3 = Enigma!3;

/// The M4 Enigma machine, which has four rotor slots. The fourth rotor never rotates.
alias EnigmaM4 = Enigma!(4, true);

/// Predefined existent rotors.
auto rotorI(dchar ringSetting = 'A') pure
{
    return rotor("EKMFLGDQVZNTOWYHXUSPAIBRCJ", 'Q', ringSetting);
}

/// ditto
auto rotorII(dchar ringSetting = 'A') pure
{
    return rotor("AJDKSIRUXBLHWTMCQGZNPYFVOE", 'E', ringSetting);
}

/// ditto
auto rotorIII(dchar ringSetting = 'A') pure
{
    return rotor("BDFHJLCPRTXVZNYEIWGAKMUSQO", 'V', ringSetting);
}

/// ditto
auto rotorIV(dchar ringSetting = 'A') pure
{
    return rotor("ESOVPZJAYQUIRHXLNFTGKDCMWB", 'J', ringSetting);
}

/// ditto
auto rotorV(dchar ringSetting = 'A') pure
{
    return rotor("VZBRGITYUPSDNHLXAWMJQOFECK", 'Z', ringSetting);
}

/// ditto
auto rotorVI(dchar ringSetting = 'A') pure
{
    return rotor("JPGVOUMFYQBENHZRDKASXLICTW", 'Z', 'M', ringSetting);
}

/// ditto
auto rotorVII(dchar ringSetting = 'A') pure
{
    return rotor("NZJHGRCXMYSWBOUFAIVLPEKQDT", 'Z', 'M', ringSetting);
}

/// ditto
auto rotorVIII(dchar ringSetting = 'A') pure
{
    return rotor("FKQHTLXOCBJSPDZRAMEWNIUYGV", 'Z', 'M', ringSetting);
}

/++
+ Predefined existent rotors. Because these rotors have no turnover notches, they are generally set
+ side by side with a reflector.
+/
auto rotorBeta(dchar ringSetting = 'A') pure
{
    return rotor("LEYJVCNIXWPBQMDRTAKZGFUHOS", ringSetting);
}

/// ditto
auto rotorGamma(dchar ringSetting = 'A') pure
{
    return rotor("FSOKANUERHMBTIYCWLQPZXVGJD", ringSetting);
}

/// Predefined the simplest plugboard which does not substitute.
auto plugboardDoNothing() pure
{
    return plugboard("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
}

/// Predefined existent reflectors.
auto reflectorB() pure
{
    return reflector("YRUHQSLDPXNGOKMIEBFZCWVJAT");
}

/// ditto
auto reflectorC() pure
{
    return reflector("FVPJIAOYEDRZXWGCTKUQSBNMHL");
}

/// ditto
auto reflectorBThin() pure
{
    return reflector("ENKQAUYWJICOPBLMDXZVFTHRGS");
}

/// ditto
auto reflectorCThin() pure
{
    return reflector("RDOBJNTKVEHMLFCWZAXGYIPSUQ");
}

// Double stepping test (http://www.cryptomuseum.com/crypto/enigma/working.htm)
unittest
{
    auto m3 = EnigmaM3(plugboardDoNothing, rotorI, rotorII, rotorIII, reflectorB, "oda");

    assert(m3.rotationStates == [14, 3, 0]);
    m3('A');
    assert(m3.rotationStates == [15, 3, 0]);
    m3('A');
    assert(m3.rotationStates == [16, 3, 0]);
    m3('A');
    assert(m3.rotationStates == [17, 4, 0]);
    m3('A');
    assert(m3.rotationStates == [18, 5, 1]);
    m3('A');
    assert(m3.rotationStates == [19, 5, 1]);
    m3('A');
    assert(m3.rotationStates == [20, 5, 1]);
}

/// Step-by-step enciphering.
unittest
{
    immutable pbCI = plugboard("ABIDEFGHCJKLMNOPQRSTUVWXYZ"); // C <-> I
    immutable refB = reflector("YRUHQSLDPXNGOKMIEBFZCWVJAT");
    immutable rot1 = rotor("EKMFLGDQVZNTOWYHXUSPAIBRCJ", 'Q', 'A');
    immutable rot2 = rotor("AJDKSIRUXBLHWTMCQGZNPYFVOE", 'E', 'B');
    immutable rot3 = rotor("BDFHJLCPRTXVZNYEIWGAKMUSQO", 'V', 'A');

    auto e3 = Enigma!3(pbCI, rot1, rot2, rot3, refB, ['X', 'Q', 'E']);
    assert(e3('A') == 'K');
    assert(e3('a') == 'T'); // A lowercase is automatically converted to an uppercase.
    assert(e3('5') == '5'); // A non-alphabetical character does not changes
    assert(e3('Ü') == 'Ü'); // the machine state and will be output as it is.
    assert(e3('A') == 'Q');
}

/// Encipherment is decipherment.
unittest
{
    // These have the same settings.
    auto encipherer = Enigma!2(plugboardDoNothing, rotorVI, rotorVII, reflectorC, "PY");
    auto decipherer = Enigma!2(plugboardDoNothing, rotorVI, rotorVII, reflectorC, "PY");

    foreach (dchar c; "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    {
        auto enciphered = encipherer(c);
        auto deciphered = decipherer(enciphered);
    }
}

/++
 + A certain equivalence of the M3 and the M4.
 +/
unittest
{
    // These have the equivalent settings.
    auto m3 = EnigmaM3(plugboardDoNothing, rotorIII, rotorII, rotorI, reflectorB /*!*/ ,
        "FOO");
    auto m4 = EnigmaM4(plugboardDoNothing, rotorIII, rotorII, rotorI,
        rotorBeta('A') /*!*/ , reflectorBThin /*!*/ , "FOOA" /*!*/ ); // FOO*A*

    // If each machine has just one movable rotor...
    auto e1 = Enigma!1(plugboardDoNothing, rotorI, reflectorC /*!*/ , "X");
    auto e2fixed = Enigma!(2, true  /*!*/ )(plugboardDoNothing, rotorI,
        rotorGamma('A') /*!*/ , reflectorCThin /*!*/ , "XA" /*!*/ ); // X*A*

    foreach (dchar c; "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    {
        assert(m3(c) == m4(c));
        assert(e1(c) == e2fixed(c));
    }
}

/// Encipher with the M4 and decipher with the equivalent M3.
unittest
{
    import std.algorithm.comparison : equal;
    import std.algorithm.iteration : each, map;
    import std.array : appender;

    // These have the equivalent settings.
    auto m4 = EnigmaM4(plugboard("SBCDEGFHIJKLMNOPQRATUVWXYZ"), rotorIII('Y'),
        rotorII('V'), rotorI('R'), rotorBeta, reflectorBThin, "UEQA");
    auto m3 = EnigmaM3(plugboard("SBCDEGFHIJKLMNOPQRATUVWXYZ"), rotorIII('Y'),
        rotorII('V'), rotorI('R'), reflectorB, "UEQ");

    auto enciphered = appender!dstring;
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ".map!m4.each!(c => enciphered.put(c));
    assert(enciphered.data == "RIIGSIBEBIZKCTZSSDGQMLSVUX");

    auto deciphered = enciphered.data.map!m3;
    assert(deciphered.equal("ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
}
