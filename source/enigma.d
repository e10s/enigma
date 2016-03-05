// Written in the D programming language.

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
    in
    {
        import boolean_matrix : isBijective;

        assert(perm.isBijective, "Rotor must be bijective.");
    }
    body
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

    /++
     + Constructs a rotor with five turnover notches.
     + If turnover1 is `1` and turnover2 is `25`, the next rotor steps
     + when this rotor steps from B to C and from Z to A.
     + If ringOffset is `2`, it corresponds to "C-03".
     +/
    this()(in auto ref BSM!N perm, size_t turnover1, size_t turnover2, size_t turnover3, size_t turnover4, size_t turnover5, size_t ringOffset) pure
    {
        import std.algorithm.sorting : sort;
        import std.array : array;

        this(perm, ringOffset);
        this.turnovers = [turnover1 % N, turnover2 % N, turnover3 % N, turnover4 % N, turnover5 % N].sort().array;
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
    assert(ringSetting.isAlpha, "Bad ring setting.");
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
    assert(turnover.isAlpha, "Bad turnover setting.");
    assert(ringSetting.isAlpha, "Bad ring setting.");
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
    assert(turnover1.isAlpha, "Bad turnover setting.");
    assert(turnover2.isAlpha, "Bad turnover setting.");
    assert(ringSetting.isAlpha, "Bad ring setting.");
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

/++
 + A convenience function to make a rotor with five turnover notches
 + from a forward substitution.
 + If turnover1 is `'B'` and turnover2 is `'Z'`, the next rotor steps
 + when this rotor steps from B to C and from Z to A.
 + If ringSetting is `'C'`, it corresponds to "C-03".
 +/
auto rotor(S)(in S forwardSubstitution, dchar turnover1, dchar turnover2, dchar turnover3, dchar turnover4, dchar turnover5, dchar ringSetting) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength;

    assert(forwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(forwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    assert(turnover1.isAlpha, "Bad turnover setting.");
    assert(turnover2.isAlpha, "Bad turnover setting.");
    assert(turnover3.isAlpha, "Bad turnover setting.");
    assert(turnover4.isAlpha, "Bad turnover setting.");
    assert(turnover5.isAlpha, "Bad turnover setting.");
    assert(ringSetting.isAlpha, "Bad ring setting.");
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Rotor(forwardSubstitution.map!toUpper.map!"a-'A'".array.permutation!N,
        turnover1.toUpper - 'A', turnover2.toUpper - 'A', turnover3.toUpper - 'A',
        turnover4.toUpper - 'A', turnover5.toUpper - 'A', ringSetting.toUpper - 'A');
}

///
struct EntryWheel
{
    import boolean_matrix : BSM;

    immutable BSM!N perm;
    ///
    this()(in auto ref BSM!N perm) pure
    in
    {
        import boolean_matrix : isBijective;

        assert(perm.isBijective, "Entry wheel must be bijective.");
    }
    body
    {
        this.perm = cast(immutable) perm;
    }

    alias perm this;
}

/++
 + A convenience function to make an entry wheel from a substitution.
 +/
auto entryWheel(S)(in S backwardSubstitution) pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : toUpper;
    import std.range : iota, walkLength;

    assert(backwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(backwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation, transpose;

    return EntryWheel(backwardSubstitution.map!toUpper.map!"a-'A'".array.permutation!N.transpose);
}

///
struct Plugboard
{
    import boolean_matrix : BSM;

    immutable BSM!N perm;
    ///
    this()(in auto ref BSM!N perm) pure
    in
    {
        import boolean_matrix : isBijective, isSymmetric;

        assert(perm.isBijective, "Plugboard must be bijective.");
        assert(perm.isSymmetric, "Plugboard must be symmetric.");
    }
    body
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
    this()(in auto ref BSM!N perm, size_t ringOffset) pure
    in
    {
        import boolean_matrix : isBijective, isIrreflexive, isSymmetric;

        assert(perm.isBijective, "Reflector must be bijective.");
        assert(perm.isIrreflexive, "Reflector must be irreflexive.");
        assert(perm.isSymmetric, "Reflector must be bijective.");
    }
    body
    {
        import boolean_matrix : lowerRotator, upperRotator;

        this.perm = lowerRotator!N(ringOffset) * perm * upperRotator!N(ringOffset);
    }

    alias perm this;
}

/++
 + A convenience function to make a reflector from a substitution.
 +/
auto reflector(S)(in S substitution, dchar ringSetting = 'A') pure if (isSomeStringOrDcharRange!S)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.algorithm.searching : canFind;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength, zip;

    assert(substitution.walkLength == N, "Bad length.");
    assert(!N.iota.zip(substitution.map!toUpper.map!"a-'A'").canFind!"a[0]==a[1]", "Self-loop found.");
    assert(N.iota.isPermutation(substitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    assert(ringSetting.isAlpha, "Bad ring setting.");
}
body
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import boolean_matrix : permutation;

    return Reflector(substitution.map!toUpper.map!"a-'A'".array.permutation!N,
        ringSetting.toUpper - 'A');
}

/// Currently machines with the double-stepping mechanism are available.
struct Enigma(size_t rotorN, bool fixedFinalRotor = false, bool hasPlugboard = true, bool settableReflectorPos = false)
{
    import boolean_matrix : BSM;
    private immutable BSM!N composedInputPerm;
    private immutable Rotor[rotorN] rotors;
    private immutable BSM!N reflector;
    private size_t[rotorN] rotationStates;

    import meta_workaround : Repeat;

    ///
    this(in EntryWheel entryWheel, in Repeat!(rotorN, Rotor) rotors,
        in Reflector reflector, in dchar[rotorN] rotorStartPos)
    in
    {
        foreach (dchar c; rotorStartPos)
        {
            import std.ascii : isAlpha;

            assert(c.isAlpha, "Bad start position.");
        }
    }
    body
    {
        foreach (i, ref e; rotationStates)
        {
            import std.ascii : toUpper;

            e = rotorStartPos[i].toUpper - 'A';
        }

        this.composedInputPerm = cast(immutable) entryWheel.perm;
        this.rotors[] = cast(immutable)[rotors][];
        this.reflector = cast(immutable) reflector.perm;
    }

    ///
    static if (settableReflectorPos)
    {
        this(in EntryWheel entryWheel, in Repeat!(rotorN, Rotor) rotors,
            in Reflector reflector, in dchar[rotorN] rotorStartPos,
            dchar reflectorPos)
        in
        {
            import std.ascii : isAlpha;

            assert(reflectorPos.isAlpha, "Bad reflector position.");
        }
        body
        {
            this(entryWheel, rotors, reflector, rotorStartPos);

            import std.ascii : toUpper;
            import boolean_matrix : lowerRotator, upperRotator;

            immutable refOffset = reflectorPos.toUpper - 'A';
            this.reflector = upperRotator!N(refOffset) * reflector * lowerRotator!N(refOffset);
        }
    }

    ///
    static if (hasPlugboard)
    {
        this(in Plugboard plugboard, in EntryWheel entryWheel,
            in Repeat!(rotorN, Rotor) rotors,
            in Reflector reflector, in dchar[rotorN] rotorStartPos)
        {
            this(entryWheel, rotors, reflector, rotorStartPos);
            this.composedInputPerm = cast(immutable) (entryWheel * plugboard);
        }

        ///
        static if (settableReflectorPos)
        {
            this(in Plugboard plugboard, in EntryWheel entryWheel,
                in Repeat!(rotorN, Rotor) rotors,
                in Reflector reflector, in dchar[rotorN] rotorStartPos,
                dchar reflectorPos)
            in
            {
                import std.ascii : isAlpha;

                assert(reflectorPos.isAlpha, "Bad reflector position.");
            }
            body
            {
                this(plugboard, entryWheel, rotors, reflector, rotorStartPos);

                import std.ascii : toUpper;
                import boolean_matrix : lowerRotator, upperRotator;

                immutable refOffset = reflectorPos.toUpper - 'A';
                this.reflector = upperRotator!N(refOffset) * reflector * lowerRotator!N(refOffset);
            }
        }
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
        import boolean_matrix : lowerRotator, upperRotator;

        immutable ptrdiff_t x = rotorID == 0 ? rotationStates[0] : rotationStates[rotorID] - rotationStates[rotorID - 1];
        immutable relRotator = x > 0 ? lowerRotator!N(x) : upperRotator!N(-x);
        immutable composedPerm = rotors[rotorID].perm * relRotator * prevPerm;
        return rotorID == rotorN - 1 ? upperRotator!N(rotationStates[rotorID]) * composedPerm
            : composeForwardPermutation(composedPerm, rotorID + 1);
    }

    private auto process(size_t keyInputID)
    out (r)
    {
        assert(r >= 0);
    }
    body
    {
        step();

        import boolean_matrix : isBijective, isIrreflexive, isSymmetric, transpose, BCV;

        immutable fwdPerm = composeForwardPermutation(composedInputPerm, 0);
        // bwdPerm = fwdPerm^-1 = fwdPerm^T
        immutable wholePerm = fwdPerm.transpose * reflector * fwdPerm;
        assert(wholePerm.isBijective);
        assert(wholePerm.isIrreflexive);
        assert(wholePerm.isSymmetric);
        immutable w = wholePerm * BCV!N.e(keyInputID);
        import std.algorithm.searching : countUntil;

        immutable r = w[].countUntil!"a";
        return r;
    }

    /// Enciphers only an alphabetical character through the current Enigma machine.
    dchar opCall(dchar keyInput)
    {
        import std.ascii : isAlpha, toUpper;

        return keyInput.isAlpha ? process(keyInput.toUpper - 'A') + 'A' : keyInput;
    }
}

/// Enigma I 'Wehrmacht', which has three rotor slots.
alias EnigmaI = Enigma!3;

/// Enigma M3, which has three rotor slots.
alias EnigmaM3 = EnigmaI;

/// Enigma M4, which has four rotor slots. The fourth rotor never rotates.
alias EnigmaM4 = Enigma!(4, true);

/// Norway Enigma, which has three rotor slots.
alias Norway = EnigmaI;

/// Enigma D, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaD = Enigma!(3, false, false, true);

/// Swiss K, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias SwissK = EnigmaD;

/// Enigma R, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaR = EnigmaD;

/// Enigma T 'Tirpitz', which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaT = EnigmaD;

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

/// ditto
auto rotorINor(dchar ringSetting = 'A') pure
{
    return rotor("WTOKASUYVRBXJHQCPZEFMDINLG", 'Q', ringSetting);
}

/// ditto
auto rotorIINor(dchar ringSetting = 'A') pure
{
    return rotor("GJLPUBSWEMCTQVHXAOFZDRKYNI", 'E', ringSetting);
}

/// ditto
auto rotorIIINor(dchar ringSetting = 'A') pure
{
    return rotor("JWFMHNBPUSDYTIXVZGRQLAOEKC", 'V', ringSetting);
}

/// ditto
alias rotorIVNor = rotorIV;

/// ditto
auto rotorVNor(dchar ringSetting = 'A') pure
{
    return rotor("HEJXQOTZBVFDASCILWPGYNMURK", 'Z', ringSetting);
}

/// ditto
auto rotorID(dchar ringSetting = 'A') pure
{
    return rotor("LPGSZMHAEOQKVXRFYBUTNICJDW", 'Y', ringSetting);
}

/// ditto
auto rotorIID(dchar ringSetting = 'A') pure
{
    return rotor("SLVGBTFXJQOHEWIRZYAMKPCNDU", 'E', ringSetting);
}

/// ditto
auto rotorIIID(dchar ringSetting = 'A') pure
{
    return rotor("CJGDPSHKTURAWZXFMYNQOBVLIE", 'N', ringSetting);
}

/// ditto
auto rotorIK(dchar ringSetting = 'A') pure
{
    return rotor("PEZUOHXSCVFMTBGLRINQJWAYDK", 'Y', ringSetting);
}

/// ditto
auto rotorIIK(dchar ringSetting = 'A') pure
{
    return rotor("ZOUESYDKFWPCIQXHMVBLGNJRAT", 'E', ringSetting);
}

/// ditto
auto rotorIIIK(dchar ringSetting = 'A') pure
{
    return rotor("EHRVXGAOBQUSIMZFLYNWKTPDJC", 'N', ringSetting);
}

/// ditto
auto rotorIR(dchar ringSetting = 'A') pure
{
    return rotor("JGDQOXUSCAMIFRVTPNEWKBLZYH", 'N', ringSetting);
}

/// ditto
auto rotorIIR(dchar ringSetting = 'A') pure
{
    return rotor("NTZPSFBOKMWRCJDIVLAEYUXHGQ", 'E', ringSetting);
}

/// ditto
auto rotorIIIR(dchar ringSetting = 'A') pure
{
    return rotor("JVIUBHTCDYAKEQZPOSGXNRMWFL", 'Y', ringSetting);
}

/// ditto
auto rotorIT(dchar ringSetting = 'A') pure
{
    return rotor("KPTYUELOCVGRFQDANJMBSWHZXI", 'W', 'Z', 'E', 'K', 'Q', ringSetting);
}

auto rotorIIT(dchar ringSetting = 'A') pure
{
    return rotor("UPHZLWEQMTDJXCAKSOIGVBYFNR", 'W', 'Z', 'F', 'L', 'R', ringSetting);
}

auto rotorIIIT(dchar ringSetting = 'A') pure
{
    return rotor("QUDLYRFEKONVZAXWHMGPJBSICT", 'W', 'Z', 'E', 'K', 'Q', ringSetting);
}

auto rotorIVT(dchar ringSetting = 'A') pure
{
    return rotor("CIWTBKXNRESPFLYDAGVHQUOJZM", 'W', 'Z', 'F', 'L', 'R', ringSetting);
}

auto rotorVT(dchar ringSetting = 'A') pure
{
    return rotor("UAXGISNJBVERDYLFZWTPCKOHMQ", 'Y', 'C', 'F', 'K', 'R', ringSetting);
}

auto rotorVIT(dchar ringSetting = 'A') pure
{
    return rotor("XFUZGALVHCNYSEWQTDMRBKPIOJ", 'X', 'E', 'I', 'M', 'Q', ringSetting);
}

auto rotorVIIT(dchar ringSetting = 'A') pure
{
    return rotor("BJVFTXPLNAYOZIKWGDQERUCHSM", 'Y', 'C', 'F', 'K', 'R', ringSetting);
}

auto rotorVIIIT(dchar ringSetting = 'A') pure
{
    return rotor("YMTPNZHWKODAJXELUQVGCBISFR", 'X', 'E', 'I', 'M', 'Q', ringSetting);
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

/// Predefined the simplest entry wheel which does not substitute.
auto entryWheelABC() pure
{
    return entryWheel("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
}

/// Predefined entry wheel: QWE... -> ABC...
auto entryWheelQWE() pure
{
    return entryWheel("QWERTZUIOASDFGHJKPYXCVBNML");
}

/// Predefined entry wheel: KZR... -> ABC...
auto entryWheelKZR() pure
{
    return entryWheel("KZROUQHYAIGBLWVSTDXFPNMCJE");
}

/// Predefined the simplest plugboard which does not substitute.
auto plugboardDoNothing() pure
{
    return plugboard("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
}

/// Predefined existent reflectors.
auto reflectorA() pure
{
    return reflector("EJMZALYXVBWFCRQUONTSPIKHGD");
}

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

/// ditto
auto reflectorNor() pure
{
    return reflector("MOWJYPUXNDSRAIBFVLKZGQCHET");
}

/// ditto
auto reflectorD(dchar ringSetting = 'A') pure
{
    return reflector("IMETCGFRAYSQBZXWLHKDVUPOJN", ringSetting);
}

/// ditto
alias reflectorK = reflectorD;

/// ditto
auto reflectorR(dchar ringSetting = 'A') pure
{
    return reflector("QYHOGNECVPUZTFDJAXWMKISRBL", ringSetting);
}

/// ditto
auto reflectorT(dchar ringSetting = 'A') pure
{
    return reflector("GEKPBTAUMOCNILJDXZYFHWVQSR", ringSetting);
}

// Double stepping test (http://www.cryptomuseum.com/crypto/enigma/working.htm)
unittest
{
    auto m3 = EnigmaM3(plugboardDoNothing, entryWheelABC, rotorI, rotorII, rotorIII, reflectorB, "ODA");

    assert(m3.rotationStates == [14, 3, 0]);

    assert(m3('A') == 'H');
    assert(m3.rotationStates == [15, 3, 0]);

    assert(m3('A') == 'D');
    assert(m3.rotationStates == [16, 3, 0]);

    assert(m3('A') == 'Z');
    assert(m3.rotationStates == [17, 4, 0]);

    assert(m3('A') == 'G');
    assert(m3.rotationStates == [18, 5, 1]);

    assert(m3('A') == 'O');
    assert(m3.rotationStates == [19, 5, 1]);

    assert(m3('A') == 'V');
    assert(m3.rotationStates == [20, 5, 1]);
}

// The noches are fixed to the ring.
unittest
{
    auto ed = EnigmaD(entryWheelQWE, rotorID, rotorIID, rotorIIID, reflectorD, "UDN" /*!*/ , 'B');

    assert(ed.rotationStates == [20, 3, 13]);

    assert(ed('A') == 'Z');
    assert(ed.rotationStates == [21, 3, 13]);

    assert(ed('A') == 'D');
    assert(ed.rotationStates == [22, 3, 13]);

    assert(ed('A') == 'V');
    assert(ed.rotationStates == [23, 3, 13]);

    assert(ed('A') == 'I');
    assert(ed.rotationStates == [24, 3, 13]);

    assert(ed('A') == 'C');
    assert(ed.rotationStates == [25, 4, 13]);

    assert(ed('A') == 'Z');
    assert(ed.rotationStates == [0, 5, 14]);


    // The K's rotor positions are same as the D's.
    auto sk = SwissK(entryWheelQWE, rotorIK('Z'), rotorIIK('Y'), rotorIIIK, reflectorK('E'), "UDN" /*!*/ , 'X');

    assert(sk.rotationStates == [20, 3, 13]);

    assert(sk('A') == 'Y');
    assert(sk.rotationStates == [21, 3, 13]);

    assert(sk('A') == 'H');
    assert(sk.rotationStates == [22, 3, 13]);

    assert(sk('A') == 'U');
    assert(sk.rotationStates == [23, 3, 13]);

    assert(sk('A') == 'M');
    assert(sk.rotationStates == [24, 3, 13]);

    assert(sk('A') == 'V');
    assert(sk.rotationStates == [25, 4, 13]);

    assert(sk('A') == 'Q');
    assert(sk.rotationStates == [0, 5, 14]);
}

unittest
{
    auto nor1 = Norway(entryWheelABC, rotorINor, rotorIINor, rotorIIINor, reflectorNor, "PDL");

    assert(nor1('A') == 'I');
    assert(nor1('A') == 'F');
    assert(nor1('A') == 'L');
    assert(nor1('A') == 'F');
    assert(nor1('A') == 'F');
    assert(nor1('A') == 'H');


    auto nor2 = Norway(plugboard("ABCDEGFHIJKLMNOPQRSTUVWXYZ"), entryWheelABC, rotorINor, rotorIVNor, rotorVNor, reflectorNor, "PDL");

    assert(nor2('A') == 'P');
    assert(nor2('A') == 'G');
    assert(nor2('A') == 'J');
    assert(nor2('A') == 'N');
    assert(nor2('A') == 'U');
    assert(nor2('A') == 'Z');
}

unittest
{
    auto er = EnigmaR(entryWheelQWE, rotorIR, rotorIIR, rotorIIIR, reflectorR('E'), "KBA", 'C');

    assert(er('A') == 'Y');
    assert(er('A') == 'P');
    assert(er('A') == 'I');
    assert(er('A') == 'R');
    assert(er('A') == 'Z');
    assert(er('A') == 'S');
}

unittest
{
    auto et1 = EnigmaT(entryWheelKZR, rotorIT, rotorIIT, rotorIIIT, reflectorT, "AFR", 'H');

    assert(et1('A') == 'E');
    assert(et1('A') == 'T');
    assert(et1('A') == 'U');
    assert(et1('A') == 'I');
    assert(et1('A') == 'H');
    assert(et1('A') == 'M');


    auto et2 = EnigmaT(entryWheelKZR, rotorVIT, rotorVT, rotorIVT, reflectorT, "AFR", 'A');

    assert(et2('A') == 'X');
    assert(et2('A') == 'F');
    assert(et2('A') == 'M');
    assert(et2('A') == 'R');
    assert(et2('A') == 'C');
    assert(et2('A') == 'B');


    auto et3 = EnigmaT(entryWheelKZR, rotorVIIT, rotorVIIIT('Z'), rotorIVT, reflectorT('C'), "AFR", 'V');

    assert(et3('A') == 'D');
    assert(et3('A') == 'S');
    assert(et3('A') == 'F');
    assert(et3('A') == 'N');
    assert(et3('A') == 'D');
    assert(et3('A') == 'T');
}

/// Step-by-step enciphering.
unittest
{
    immutable pbCI = plugboard("ABIDEFGHCJKLMNOPQRSTUVWXYZ"); // C <-> I
    immutable enWh = entryWheel("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    immutable refB = reflector("YRUHQSLDPXNGOKMIEBFZCWVJAT");
    immutable rot1 = rotor("EKMFLGDQVZNTOWYHXUSPAIBRCJ", 'Q', 'A');
    immutable rot2 = rotor("AJDKSIRUXBLHWTMCQGZNPYFVOE", 'E', 'B');
    immutable rot3 = rotor("BDFHJLCPRTXVZNYEIWGAKMUSQO", 'V', 'A');

    auto e3 = Enigma!3(pbCI, enWh, rot1, rot2, rot3, refB, ['X', 'Q', 'E']);
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
    auto encipherer = Enigma!2(entryWheelABC, rotorVI, rotorVII, reflectorC, "PY");
    auto decipherer = Enigma!2(entryWheelABC, rotorVI, rotorVII, reflectorC, "PY");

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
    auto m3 = EnigmaM3(entryWheelABC, rotorIII, rotorII, rotorI, reflectorB /*!*/ ,
        "FOO");
    auto m4 = EnigmaM4(entryWheelABC, rotorIII, rotorII, rotorI,
        rotorBeta('A') /*!*/ , reflectorBThin /*!*/ , "FOOA" /*!*/ ); // FOO*A*

    // If each machine has just one movable rotor...
    auto e1 = Enigma!1(entryWheelABC, rotorI, reflectorC /*!*/ , "X");
    auto e2fixed = Enigma!(2, true  /*!*/ )(entryWheelABC, rotorI,
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
    auto m4 = EnigmaM4(plugboard("SBCDEGFHIJKLMNOPQRATUVWXYZ"), entryWheelABC, rotorIII('Y'),
        rotorII('V'), rotorI('R'), rotorBeta, reflectorBThin, "UEQA");
    auto m3 = EnigmaM3(plugboard("SBCDEGFHIJKLMNOPQRATUVWXYZ"), entryWheelABC, rotorIII('Y'),
        rotorII('V'), rotorI('R'), reflectorB, "UEQ");

    auto enciphered = appender!dstring;
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ".map!m4.each!(c => enciphered.put(c));
    assert(enciphered.data == "RIIGSIBEBIZKCTZSSDGQMLSVUX");

    auto deciphered = enciphered.data.map!m3;
    assert(deciphered.equal("ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
}
