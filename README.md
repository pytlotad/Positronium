# Klasyczny atom wodoru — symulacja elektrodynamiczna

Program w CERN ROOT przedstawia klasyczny model elektrodynamiczny układu
elektron–proton:

- oba ładunki poruszają się w polu Coulomba w układzie środka masy;
- całkowanie ruchu wykorzystuje metodę velocity-Verlet;
- promieniowanie Larmora jest liczone dla obu przyspieszających ładunków;
- energia wypromieniowana jest odejmowana od ruchu względnego, więc elektron
  spiralnie zbliża się do protonu.

Każda cząstka ma też klasyczny moment dipolowy o losowym (lecz powtarzalnym)
kierunku początkowym. Symulacja uwzględnia pole dipola, siłę dipol–dipol oraz
moment skręcający, który zmienia kierunki strzałek. Rzeczywista precesja jest
na tej skali czasu zbyt wolna do zauważenia, dlatego wyłącznie wizualna ewolucja
orientacji dipoli jest przyspieszona; ich wpływ mechaniczny nie jest wzmacniany.

To jest oczekiwany rezultat klasycznej elektrodynamiki: przyspieszany ładunek
promieniuje, dlatego taka orbita traci energię. Stabilność prawdziwego wodoru
opisuje mechanika kwantowa. Program kończy ruch przy `1e-14 m`, ponieważ model
punktowych ładunków nie jest już tam wiarygodny.

## Wymagania i uruchomienie

Wymagany jest CERN ROOT (z `root-config` dostępnym w `PATH`).

```bash
make
./atom_sim
```

Wykres jest w jednostkach promienia Bohra `a₀`. Odczyt na dole zawiera czas w
pikosekundach oraz całkowitą energię mechaniczną układu elektron–proton w eV.
Wyświetlana jest również skumulowana energia wypromieniowana w eV. Linia nad
bieżącymi licznikami pokazuje ich wartości początkowe, a linia pod nimi zmiany
względem początku, oznaczone symbolem `Δ`.

Przyciski nad widokiem: `STOP` wstrzymuje animację (zmienia się na `START`), a
`EXIT` natychmiast zamyka program.

Symulacja jest renderowana w 3D. Przeciągnij myszą w obszarze orbity, aby
zmienić kąt kamery i obejrzeć ruch z innej perspektywy.
