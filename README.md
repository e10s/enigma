# enigma - Enigma machine simulation library

enigma is a library for simulating the Enigma machines written in D. It aims to reproduce the existent variations of machines and mechanisms and be implemented using the mathematical representation of the Enigma.

https://en.wikipedia.org/wiki/Enigma_machine

## Example

```d
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
```
