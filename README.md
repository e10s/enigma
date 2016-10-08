# enigma - Enigma machine simulation library

*enigma* is a library for simulating the Enigma machines written in the D programming language.

This library is aims to reproduce the existent variations of machines and mechanisms, and be also able to create user-defined models.
Moreover, it is intended to be implemented internally using the mathematical representation of the components, such as groups and matrices.

**Note that of course such replicas are extremely vulnerable for actual ciphering use!**

## Example

This simple program using *enigma* imitates the existent "M3" model. It will get user input and write the ciphered text to `stdout`.

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

Most of the known models and their components are predefined as seen above.
For more information, read the documents generated with `dub --build=docs`.

## About Enigma
* [Enigma machine on Wikipedia](https://en.wikipedia.org/wiki/Enigma_machine)
* [Enigma Cipher Machines on Crypto Museum](http://www.cryptomuseum.com/crypto/enigma/index.htm)
