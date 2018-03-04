// Written in the D programming language.

/*
Copyright Kazuya Takahashi 2016-2018.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
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
    import structure : PermutationElement;

    immutable PermutationElement!N perm;
    private immutable bool hasNotch = false;
    private immutable size_t[] turnovers;

    import std.meta : allSatisfy;
    import std.traits : isIntegral;
    /++
     + Constructs a rotor.
     + If turnovers are not specified, the rotor has no turnover notches,
     + otherwise, for example, if turnovers are `1` and `25`, the next rotor
     + steps when this rotor steps from B to C and from Z to A.
     + If ringOffset is `2`, it corresponds to "C-03".
     +/
    this(I...)(in auto ref PermutationElement!N perm, I turnovers, size_t ringOffset) pure
        if (allSatisfy!(isIntegral, I) && I.length <= N)
    in
    {
        assert(perm.isBijective, "Rotor must be bijective.");
        foreach (t; turnovers)
        {
            assert(t >= 0, "Turnover must be positive.");
        }
    }
    do
    {
        import std.algorithm.iteration : map, uniq;
        import std.algorithm.sorting : sort;
        import std.array : array;
        import structure : cyclicPermutation, cyclicPermutationInv;

        this.perm = cyclicPermutationInv!N(ringOffset) * perm * cyclicPermutation!N(ringOffset);
        size_t[] ts = [turnovers];
        this.turnovers = ts.map!(a => a % N).array.sort().uniq.array.idup;
        hasNotch = I.length > 0;
    }
}

import std.meta : allSatisfy;
import std.traits : isSomeChar;
/++
 + A convenience function to make a rotor from a forward substitution.
 + If turnovers are not specified, the rotor has no turnover notches,
 + otherwise, for example, if  turnovers are `'B'` and `'Z'`, the next rotor
 + steps when this rotor steps from B to C and from Z to A.
 + If ringSetting is `'C'`, it corresponds to "C-03".
 +/
auto rotor(S, C...)(in S forwardSubstitution, C turnovers, dchar ringSetting) pure
    if (isSomeStringOrDcharRange!S && allSatisfy!(isSomeChar, C) && C.length <= N)
in
{
    import std.algorithm.comparison : isPermutation;
    import std.algorithm.iteration : map;
    import std.ascii : isAlpha, toUpper;
    import std.range : iota, walkLength;

    assert(forwardSubstitution.walkLength == N, "Bad length.");
    assert(N.iota.isPermutation(forwardSubstitution.map!toUpper.map!"a-'A'"), "Bad permutation.");
    foreach (dchar t; turnovers)
    {
        assert(t.isAlpha, "Bad turnover setting.");
    }
    assert(ringSetting.isAlpha, "Bad ring setting.");
}
do
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import structure : permutation;

    static if (C.length)
    {
        import std.meta : Repeat;
        import std.typecons : Tuple;

        Tuple!(Repeat!(C.length, size_t)) ts;
        foreach (i, dchar t; turnovers)
        {
            ts[i] = t.toUpper - 'A';
        }
        return Rotor(forwardSubstitution.map!toUpper.map!"a-size_t('A')".array.permutation!N,
            ts.expand, ringSetting.toUpper - 'A');
    }
    else
    {
        return Rotor(forwardSubstitution.map!toUpper.map!"a-size_t('A')".array.permutation!N,
            ringSetting.toUpper - 'A');
    }
}

///
struct EntryWheel
{
    import structure : PermutationElement;

    immutable PermutationElement!N perm;
    ///
    this()(in auto ref PermutationElement!N perm) pure
    in
    {
        assert(perm.isBijective, "Entry wheel must be bijective.");
    }
    do
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
do
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import structure : permutation;

    return EntryWheel(backwardSubstitution.map!toUpper.map!"a-size_t('A')".array.permutation!N.transpose);
}

///
struct Plugboard
{
    import structure : PermutationElement;

    immutable PermutationElement!N perm;
    ///
    this()(in auto ref PermutationElement!N perm) pure
    in
    {
        assert(perm.isBijective, "Plugboard must be bijective.");
        assert(perm.isSymmetric, "Plugboard must be symmetric.");
    }
    do
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
do
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import structure : permutation;

    return Plugboard(substitution.map!toUpper.map!"a-size_t('A')".array.permutation!N);
}

///
struct Reflector
{
    import structure : PermutationElement;

    immutable PermutationElement!N perm;
    ///
    this()(in auto ref PermutationElement!N perm, size_t ringOffset) pure
    in
    {
        assert(perm.isBijective, "Reflector must be bijective.");
        assert(perm.isIrreflexive, "Reflector must be irreflexive.");
        assert(perm.isSymmetric, "Reflector must be bijective.");
    }
    do
    {
        import structure : cyclicPermutation, cyclicPermutationInv;

        this.perm = cyclicPermutationInv!N(ringOffset) * perm * cyclicPermutation!N(ringOffset);
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
do
{
    import std.algorithm.iteration : map;
    import std.array : array;
    import std.ascii : toUpper;
    import structure : permutation;

    return Reflector(substitution.map!toUpper.map!"a-size_t('A')".array.permutation!N,
        ringSetting.toUpper - 'A');
}

///
enum EnigmaType : uint
{
    ///
    none,
    /// No double stepping.
    gearDrive = 1 << 0,
    ///
    plugboard = 1 << 1,
    ///
    fixedFinalRotor = 1 << 2,
    ///
    settableReflectorPos = 1 << 3,
    ///
    movableReflector = 1 << 4
}

///
struct Enigma(size_t rotorN, EnigmaType enigmaType = EnigmaType.none)
{
    import structure : PermutationElement;
    private immutable PermutationElement!N composedInputPerm;
    private immutable Rotor[rotorN] rotors;
    private immutable PermutationElement!N reflector;
    private size_t[rotorN + 1] rotationStates;

    import std.meta : Repeat;

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
    do
    {
        foreach (i, ref e; rotorStartPos)
        {
            import std.ascii : toUpper;

            rotationStates[i] = e.toUpper - 'A';
        }

        this.composedInputPerm = cast(immutable) entryWheel.perm;
        this.rotors[] = cast(immutable)[rotors][];
        this.reflector = cast(immutable) reflector.perm;
    }

    ///
    this()(in EntryWheel entryWheel, in Repeat!(rotorN, Rotor) rotors,
            in Reflector reflector, in dchar[rotorN] rotorStartPos, dchar reflectorPos)
            if (enigmaType & EnigmaType.settableReflectorPos
                || enigmaType & EnigmaType.movableReflector)
    in
    {
        import std.ascii : isAlpha;

        assert(reflectorPos.isAlpha, "Bad reflector position.");
    }
    do
    {
        this(entryWheel, rotors, reflector, rotorStartPos);

        import std.ascii : toUpper;

        rotationStates[$ - 1] = reflectorPos.toUpper - 'A';
    }

    ///
    this()(in Plugboard plugboard, in EntryWheel entryWheel, in Repeat!(rotorN,
            Rotor) rotors, in Reflector reflector, in dchar[rotorN] rotorStartPos)
            if (enigmaType & EnigmaType.plugboard)
    {
        this(entryWheel, rotors, reflector, rotorStartPos);
        this.composedInputPerm = cast(immutable)(entryWheel * plugboard);
    }

    ///
    this()(in Plugboard plugboard, in EntryWheel entryWheel, in Repeat!(rotorN, Rotor) rotors,
            in Reflector reflector, in dchar[rotorN] rotorStartPos, dchar reflectorPos)
            if (enigmaType & EnigmaType.plugboard && (enigmaType & EnigmaType.settableReflectorPos
                || enigmaType & EnigmaType.movableReflector))
    in
    {
        import std.ascii : isAlpha;

        assert(reflectorPos.isAlpha, "Bad reflector position.");
    }
    do
    {
        this(plugboard, entryWheel, rotors, reflector, rotorStartPos);

        import std.ascii : toUpper;

        rotationStates[$ - 1] = reflectorPos.toUpper - 'A';
    }

    private void step()
    {
        enum movableRotorN = (enigmaType & EnigmaType.fixedFinalRotor) ? rotorN - 1 :
            (enigmaType & EnigmaType.movableReflector) ? rotorN + 1 : rotorN;
        bool[movableRotorN] stepFlag;

        stepFlag[0] = true;

        foreach (rotorID; 0 .. movableRotorN - 1)
        {
            import std.algorithm.searching : canFind;

            if (rotors[rotorID].turnovers.canFind(rotationStates[rotorID]))
            {
                stepFlag[rotorID] = true;
                stepFlag[rotorID + 1] = true;
            }
            else static if (enigmaType & EnigmaType.gearDrive)
            {
                break;
            }
        }

        foreach (rotorID, e; stepFlag)
        {
            if (e)
            {
                rotationStates[rotorID] = (rotationStates[rotorID] + 1) % N;
            }
            else static if (enigmaType & EnigmaType.gearDrive)
            {
                break;
            }
        }
    }

    import structure : SetElement;

    /+
     + fwdPerm = (Um*Rm*Lm)*(Um-1*Rm-1*Lm-1)*...*(U0*R0*L0)*P
     +         = Um*(Rm*Lm*Um-1)*(Rm-1*Lm-1*Um-2)*...*(R0*L0)*P
     +         = Um*(Rm*RELm)*(Rm-1*RELm-1)*...*(R0*L0)*P
     +/
    private SetElement!N composeForward(in ref SetElement!N inputVec, size_t rotorID)
    {
        import structure : cyclicPermutation, cyclicPermutationInv;

        immutable ptrdiff_t x = rotorID == 0 ? rotationStates[0] : rotationStates[rotorID] - rotationStates[rotorID - 1];
        immutable relRotator = x > 0 ? cyclicPermutationInv!N(x) : cyclicPermutation!N(-x);
        immutable composedVec = rotors[rotorID].perm * (relRotator * inputVec);
        return rotorID == rotorN - 1 ? cyclicPermutation!N(rotationStates[rotorID]) * composedVec
            : composeForward(composedVec, rotorID + 1);
    }

    /+
     + bwdPerm = [(Um*Rm*Lm)*(Um-1*Rm-1*Lm-1)*...*(U0*R0*L0)*P]^-1
     +         = [Um*(Rm*RELm)*(Rm-1*RELm-1)*...*(R0*L0)*P]^-1
     +         = P^-1*(R0*L0)^-1*...(Rm-1*RELm-1)^-1*(Rm*RELm)^-1*Um^-1
     +         = P^-1*(U0*R0^-1)*...(RELm-1^-1*Rm-1^-1)*(RELm^-1*Rm^-1)*Lm
     +/
    private SetElement!N composeBackward(in ref SetElement!N inputVec, size_t rotorID)
    {
        import structure : cyclicPermutation, cyclicPermutationInv;

        immutable ptrdiff_t x = rotorID == 0 ? rotationStates[0] : rotationStates[rotorID] - rotationStates[rotorID - 1];
        immutable relRotatorInv = x > 0 ? cyclicPermutation!N(x) : cyclicPermutationInv!N(-x);
        immutable iv = rotorID == rotorN - 1 ? cyclicPermutationInv!N(rotationStates[rotorN - 1]) * inputVec : inputVec;
        immutable composedVec = relRotatorInv * (rotors[rotorID].perm.transpose * iv);
        return rotorID == 0 ? composedVec : composeBackward(composedVec, rotorID - 1);
    }

    // Unless the reflector is movable, the return value is constant.
    private auto composedReflector() @property
    {
        import structure : cyclicPermutation, cyclicPermutationInv;

        return cyclicPermutation!N(rotationStates[$ - 1]) * reflector * cyclicPermutationInv!N(rotationStates[$ - 1]);
    }

    private auto process(size_t keyInputID)
    out (r)
    {
        assert(r >= 0);
    }
    do
    {
        step();

        import structure : SetElement;

        immutable composedInput = composedInputPerm * SetElement!N(keyInputID);
        immutable reflectorOutput = composedReflector * composeForward(composedInput, 0);
        immutable composedOutput = composedInputPerm.transpose * composeBackward(reflectorOutput,rotorN - 1);

        return composedOutput.id;
    }

    /// Enciphers only an alphabetical character through the current Enigma machine.
    dchar opCall(dchar keyInput)
    {
        import std.ascii : isAlpha, toUpper;

        return keyInput.isAlpha ? cast(dchar)process(keyInput.toUpper - 'A') + 'A' : keyInput;
    }
}

/// Enigma I 'Wehrmacht', which has three rotor slots.
alias EnigmaI = Enigma!(3, EnigmaType.plugboard);

/// Enigma M3, which has three rotor slots.
alias EnigmaM3 = EnigmaI;

/// Enigma M4, which has four rotor slots. The fourth rotor never rotates.
alias EnigmaM4 = Enigma!(4, EnigmaType.fixedFinalRotor | EnigmaType.plugboard);

/// Norway Enigma, which has three rotor slots.
alias Norway = EnigmaI;

/// Enigma D, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaD = Enigma!(3, EnigmaType.settableReflectorPos);

/// Swiss K, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias SwissK = EnigmaD;

/// Enigma R, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaR = EnigmaD;

/// Enigma T 'Tirpitz', which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaT = EnigmaD;

/// Enigma KD, which has three rotor slots and no plugboard. The reflector can be set to any positions.
alias EnigmaKD = EnigmaD;

/// Enigma A28, which has three rotor slots and no plugboard. The reflector can be set to any positions. In addition, the rotors and the reflector rotates normally (not "double stepping").
alias EnigmaA28 = Enigma!(3, EnigmaType.gearDrive | EnigmaType.movableReflector);

/// Enigma G, which has three rotor slots and no plugboard. The reflector can be set to any positions. In addition, the rotors and the reflector rotates normally (not "double stepping").
alias EnigmaG = EnigmaA28;

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

/// ditto
auto rotorIIT(dchar ringSetting = 'A') pure
{
    return rotor("UPHZLWEQMTDJXCAKSOIGVBYFNR", 'W', 'Z', 'F', 'L', 'R', ringSetting);
}

/// ditto
auto rotorIIIT(dchar ringSetting = 'A') pure
{
    return rotor("QUDLYRFEKONVZAXWHMGPJBSICT", 'W', 'Z', 'E', 'K', 'Q', ringSetting);
}

/// ditto
auto rotorIVT(dchar ringSetting = 'A') pure
{
    return rotor("CIWTBKXNRESPFLYDAGVHQUOJZM", 'W', 'Z', 'F', 'L', 'R', ringSetting);
}

/// ditto
auto rotorVT(dchar ringSetting = 'A') pure
{
    return rotor("UAXGISNJBVERDYLFZWTPCKOHMQ", 'Y', 'C', 'F', 'K', 'R', ringSetting);
}

/// ditto
auto rotorVIT(dchar ringSetting = 'A') pure
{
    return rotor("XFUZGALVHCNYSEWQTDMRBKPIOJ", 'X', 'E', 'I', 'M', 'Q', ringSetting);
}

/// ditto
auto rotorVIIT(dchar ringSetting = 'A') pure
{
    return rotor("BJVFTXPLNAYOZIKWGDQERUCHSM", 'Y', 'C', 'F', 'K', 'R', ringSetting);
}

/// ditto
auto rotorVIIIT(dchar ringSetting = 'A') pure
{
    return rotor("YMTPNZHWKODAJXELUQVGCBISFR", 'X', 'E', 'I', 'M', 'Q', ringSetting);
}

/// ditto
auto rotorIKD(dchar ringSetting = 'A') pure
{
    return rotor("VEZIOJCXKYDUNTWAPLQGBHSFMR", 'S', 'U', 'Y', 'A', 'E', 'H', 'L', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIIKD(dchar ringSetting = 'A') pure
{
    return rotor("HGRBSJZETDLVPMQYCXAOKINFUW", 'S', 'U', 'Y', 'A', 'E', 'H', 'L', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIIIKD(dchar ringSetting = 'A') pure
{
    return rotor("NWLHXGRBYOJSAZDVTPKFQMEUIC", 'S', 'U', 'Y', 'A', 'E', 'H', 'L', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIA865(dchar ringSetting = 'A') pure
{
    return rotor("LPGSZMHAEOQKVXRFYBUTNICJDW", 'S', 'U', 'V', 'W', 'Z', 'A', 'B', 'C', 'E', 'F', 'G', 'I', 'K', 'L', 'O', 'P', 'Q', ringSetting);
}

/// ditto
auto rotorIIA865(dchar ringSetting = 'A') pure
{
    return rotor("SLVGBTFXJQOHEWIRZYAMKPCNDU", 'S', 'T', 'V', 'Y', 'Z', 'A', 'C', 'D', 'F', 'G', 'H', 'K', 'M', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIIIA865(dchar ringSetting = 'A') pure
{
    return rotor("CJGDPSHKTURAWZXFMYNQOBVLIE", 'U', 'W', 'X', 'A', 'E', 'F', 'H', 'K', 'M', 'N', 'R', ringSetting);
}

/// ditto
auto rotorIG260(dchar ringSetting = 'A') pure
{
    return rotor("RCSPBLKQAUMHWYTIFZVGOJNEXD", 'S', 'U', 'V', 'W', 'Z', 'A', 'B', 'C', 'E', 'F', 'G', 'I', 'K', 'L', 'O', 'P', 'Q', ringSetting);
}

/// ditto
auto rotorIIG260(dchar ringSetting = 'A') pure
{
    return rotor("WCMIBVPJXAROSGNDLZKEYHUFQT", 'S', 'T', 'V', 'Y', 'Z', 'A', 'C', 'D', 'F', 'G', 'H', 'K', 'M', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIIIG260(dchar ringSetting = 'A') pure
{
    return rotor("FVDHZELSQMAXOKYIWPGCBUJTNR", 'U', 'W', 'X', 'A', 'E', 'F', 'H', 'K', 'M', 'N', 'R', ringSetting);
}

/// ditto
auto rotorIG312(dchar ringSetting = 'A') pure
{
    return rotor("DMTWSILRUYQNKFEJCAZBPGXOHV", 'S', 'U', 'V', 'W', 'Z', 'A', 'B', 'C', 'E', 'F', 'G', 'I', 'K', 'L', 'O', 'P', 'Q', ringSetting);
}

/// ditto
auto rotorIIG312(dchar ringSetting = 'A') pure
{
    return rotor("HQZGPJTMOBLNCIFDYAWVEUSRKX", 'S', 'T', 'V', 'Y', 'Z', 'A', 'C', 'D', 'F', 'G', 'H', 'K', 'M', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorIIIG312(dchar ringSetting = 'A') pure
{
    return rotor("UQNTLSZFMREHDPXKIBVYGJCWOA", 'U', 'W', 'X', 'A', 'E', 'F', 'H', 'K', 'M', 'N', 'R', ringSetting);
}

/// ditto
auto rotorIG111(dchar ringSetting = 'A') pure
{
    return rotor("WLRHBQUNDKJCZSEXOTMAGYFPVI", 'S', 'U', 'V', 'W', 'Z', 'A', 'B', 'C', 'E', 'F', 'G', 'I', 'K', 'L', 'O', 'P', 'Q', ringSetting);
}

/// ditto
auto rotorIIG111(dchar ringSetting = 'A') pure
{
    return rotor("TFJQAZWMHLCUIXRDYGOEVBNSKP", 'S', 'T', 'V', 'Y', 'Z', 'A', 'C', 'D', 'F', 'G', 'H', 'K', 'M', 'N', 'Q', ringSetting);
}

/// ditto
auto rotorVG111(dchar ringSetting = 'A') pure
{
    return rotor("QTPIXWVDFRMUSLJOHCANEZKYBG", 'S', 'W', 'Z', 'F', 'H', 'M', 'Q', ringSetting);
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

/// ditto
auto reflectorFRA(dchar ringSetting = 'A') pure
{
    return reflector("KOTVPNLMJIAGHFBEWYXCZDQSRU", ringSetting);
}

/// ditto
alias reflectorA865 = reflectorD;

/// ditto
alias reflectorG260 = reflectorD;

/// ditto
auto reflectorG312(dchar ringSetting = 'A') pure
{
    return reflector("RULQMZJSYGOCETKWDAHNBXPVIF", ringSetting);
}

/// ditto
alias reflectorG111 = reflectorD;

// Double stepping test (http://www.cryptomuseum.com/crypto/enigma/working.htm)
unittest
{
    auto m3 = EnigmaM3(plugboardDoNothing, entryWheelABC, rotorI, rotorII, rotorIII, reflectorB, "ODA");

    assert(m3.rotationStates == [14, 3, 0, 0]);

    assert(m3('A') == 'H');
    assert(m3.rotationStates == [15, 3, 0, 0]);

    assert(m3('A') == 'D');
    assert(m3.rotationStates == [16, 3, 0, 0]);

    assert(m3('A') == 'Z');
    assert(m3.rotationStates == [17, 4, 0, 0]);

    assert(m3('A') == 'G');
    assert(m3.rotationStates == [18, 5, 1, 0]);

    assert(m3('A') == 'O');
    assert(m3.rotationStates == [19, 5, 1, 0]);

    assert(m3('A') == 'V');
    assert(m3.rotationStates == [20, 5, 1, 0]);
}

// The noches are fixed to the ring.
unittest
{
    auto ed = EnigmaD(entryWheelQWE, rotorID, rotorIID, rotorIIID, reflectorD, "UDN" /*!*/ , 'B');

    assert(ed.rotationStates == [20, 3, 13, 1]);

    assert(ed('A') == 'Z');
    assert(ed.rotationStates == [21, 3, 13, 1]);

    assert(ed('A') == 'D');
    assert(ed.rotationStates == [22, 3, 13, 1]);

    assert(ed('A') == 'V');
    assert(ed.rotationStates == [23, 3, 13, 1]);

    assert(ed('A') == 'I');
    assert(ed.rotationStates == [24, 3, 13, 1]);

    assert(ed('A') == 'C');
    assert(ed.rotationStates == [25, 4, 13, 1]);

    assert(ed('A') == 'Z');
    assert(ed.rotationStates == [0, 5, 14, 1]);


    // The K's rotor positions are same as the D's.
    auto sk = SwissK(entryWheelQWE, rotorIK('Z'), rotorIIK('Y'), rotorIIIK, reflectorK('E'), "UDN" /*!*/ , 'X');

    assert(sk.rotationStates == [20, 3, 13, 23]);

    assert(sk('A') == 'Y');
    assert(sk.rotationStates == [21, 3, 13, 23]);

    assert(sk('A') == 'H');
    assert(sk.rotationStates == [22, 3, 13, 23]);

    assert(sk('A') == 'U');
    assert(sk.rotationStates == [23, 3, 13, 23]);

    assert(sk('A') == 'M');
    assert(sk.rotationStates == [24, 3, 13, 23]);

    assert(sk('A') == 'V');
    assert(sk.rotationStates == [25, 4, 13, 23]);

    assert(sk('A') == 'Q');
    assert(sk.rotationStates == [0, 5, 14, 23]);
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

unittest
{
    auto ekd = EnigmaKD(entryWheelQWE, rotorIKD, rotorIIKD, rotorIIIKD, reflectorFRA, "AAA");

    assert(ekd('A') == 'W');
    assert(ekd('A') == 'C');
    assert(ekd('A') == 'H');
    assert(ekd('A') == 'U');
    assert(ekd('A') == 'G');
    assert(ekd('A') == 'K');
}

// Normal stepping and movable reflector test
unittest
{
    auto ea865 = EnigmaA28(entryWheelQWE, rotorIIIA865, rotorIIA865, rotorIA865, reflectorA865, "CAA", 'A');

    assert(ea865.rotationStates == [2, 0, 0, 0]);

    assert(ea865('A') == 'Q');
    assert(ea865.rotationStates == [3, 0, 0, 0]);

    assert(ea865('A') == 'I');
    assert(ea865.rotationStates == [4, 0, 0, 0]);

    assert(ea865('A') == 'C');
    assert(ea865.rotationStates == [5, 1, 1, 1]);

    assert(ea865('A') == 'Y');
    assert(ea865.rotationStates == [6, 2, 1, 1]);

    assert(ea865('A') == 'D');
    assert(ea865.rotationStates == [7, 2, 1, 1]);

    assert(ea865('A') == 'E');
    assert(ea865.rotationStates == [8, 3, 2, 2]);

    assert(ea865('A') == 'D');
    assert(ea865.rotationStates == [9, 3, 2, 2]);

    assert(ea865('A') == 'Z');
    assert(ea865.rotationStates == [10, 3, 2, 2]);

    assert(ea865('A') == 'X');
    assert(ea865.rotationStates == [11, 4, 3, 3]);
}

unittest
{
    auto eg260 = EnigmaG(entryWheelQWE, rotorIG260, rotorIIG260, rotorIIIG260, reflectorG260, "AAA", 'A');

    assert(eg260.rotationStates == [0, 0, 0, 0]);

    assert(eg260('A') == 'C');
    assert(eg260.rotationStates == [1, 1, 1, 1]);

    assert(eg260('A') == 'K');
    assert(eg260.rotationStates == [2, 2, 1, 1]);

    assert(eg260('A') == 'N');
    assert(eg260.rotationStates == [3, 3, 2, 1]);

    assert(eg260('A') == 'U');
    assert(eg260.rotationStates == [4, 3, 2, 1]);

    assert(eg260('A') == 'L');
    assert(eg260.rotationStates == [5, 4, 3, 1]);

    assert(eg260('A') == 'Q');
    assert(eg260.rotationStates == [6, 5, 3, 1]);

    assert(eg260('A') == 'Y');
    assert(eg260.rotationStates == [7, 6, 4, 1]);

    assert(eg260('A') == 'P');
    assert(eg260.rotationStates == [8, 6, 4, 1]);

    assert(eg260('A') == 'I');
    assert(eg260.rotationStates == [9, 7, 5, 2]);
}

unittest
{
    auto eg312 = EnigmaG(entryWheelQWE, rotorIG312, rotorIIG312('B'), rotorIIIG312, reflectorG312('Y'), "CZB", 'D');

    assert(eg312.rotationStates == [2, 25, 1, 3]);

    assert(eg312('A') == 'J');
    assert(eg312.rotationStates == [3, 0, 2, 3]);

    assert(eg312('A') == 'B');
    assert(eg312.rotationStates == [4, 0, 2, 3]);

    assert(eg312('A') == 'X');
    assert(eg312.rotationStates == [5, 1, 3, 3]);

    assert(eg312('A') == 'K');
    assert(eg312.rotationStates == [6, 2, 3, 3]);

    assert(eg312('A') == 'N');
    assert(eg312.rotationStates == [7, 3, 4, 3]);

    assert(eg312('A') == 'N');
    assert(eg312.rotationStates == [8, 3, 4, 3]);

    assert(eg312('A') == 'C');
    assert(eg312.rotationStates == [9, 4, 5, 4]);
}

unittest
{
    auto eg111 = EnigmaG(entryWheelQWE, rotorIG111, rotorVG111, rotorIIG111, reflectorG111, "AAA", 'A');

    assert(eg111.rotationStates == [0, 0, 0, 0]);

    assert(eg111('A') == 'E');
    assert(eg111.rotationStates == [1, 1, 0, 0]);

    assert(eg111('A') == 'L');
    assert(eg111.rotationStates == [2, 2, 0, 0]);

    assert(eg111('A') == 'U');
    assert(eg111.rotationStates == [3, 3, 0, 0]);

    assert(eg111('A') == 'X');
    assert(eg111.rotationStates == [4, 3, 0, 0]);

    assert(eg111('A') == 'L');
    assert(eg111.rotationStates == [5, 4, 0, 0]);

    assert(eg111('A') == 'U');
    assert(eg111.rotationStates == [6, 5, 0, 0]);

    assert(eg111('A') == 'P');
    assert(eg111.rotationStates == [7, 6, 1, 1]);

    assert(eg111('A') == 'E');
    assert(eg111.rotationStates == [8, 6, 1, 1]);

    assert(eg111('A') == 'Q');
    assert(eg111.rotationStates == [9, 7, 1, 1]);
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

    auto e3 = Enigma!(3, EnigmaType.plugboard)(pbCI, enWh, rot1, rot2, rot3, refB, ['X', 'Q', 'E']);
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
    auto e2fixed = Enigma!(2, EnigmaType.fixedFinalRotor  /*!*/ )(entryWheelABC, rotorI,
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
