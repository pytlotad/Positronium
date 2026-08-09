# Klasyczne pozytonium — symulacja elektrodynamiczna

Program w CERN ROOT przedstawia klasyczny model elektrodynamiczny układu
elektron–pozyton:

- wzajemne pola ładunków są liczone z opóźnionych pól Liénarda–Wiecherta;
- ruch jest całkowany dla relatywistycznego pędu `p = γmv`;
- reakcja promieniowania jest lokalną, zredukowaną aproksymacją
  Landaua–Lifshitza, wykorzystującą zmianę siły zewnętrznej między krokami;
- skumulowana energia promieniowania jest liczona wzorem Liénarda.

Każda cząstka ma też klasyczny moment dipolowy o losowym (lecz powtarzalnym)
kierunku początkowym. Symulacja uwzględnia pole dipola, siłę dipol–dipol oraz
moment skręcający, który zmienia kierunki strzałek. Rzeczywista precesja jest
na tej skali czasu zbyt wolna do zauważenia, dlatego wyłącznie wizualna ewolucja
orientacji dipoli jest przyspieszona; ich wpływ mechaniczny nie jest wzmacniany.

To nadal jest przybliżony model punktowych ładunków: dipole reprezentują
fenomenologicznie spin, a lokalna postać reakcji promieniowania nie zastępuje
pełnego problemu samopola. Stabilność i anihilację prawdziwego pozytonium
opisuje mechanika kwantowa. Program kończy ruch przy `1e-14 m`, ponieważ model
punktowych ładunków nie jest już tam wiarygodny.

## Wymagania i uruchomienie

Wymagany jest CERN ROOT (z `root-config` dostępnym w `PATH`).

```bash
make
```

Wykres jest w jednostkach promienia Bohra `a₀`. Odczyt na dole zawiera czas w
pikosekundach, odległość cząstek w pm oraz całkowitą energię mechaniczną układu
elektron–pozyton w eV.
Wyświetlana jest również skumulowana energia wypromieniowana w eV. Linia nad
bieżącymi licznikami pokazuje ich wartości początkowe, a linia pod nimi zmiany
względem początku, oznaczone symbolem `Δ`.

Przyciski nad widokiem: `STOP` wstrzymuje animację (zmienia się na `START`), a
`EXIT` natychmiast zamyka program.

Symulacja jest renderowana w 3D. Przeciągnij myszą w obszarze orbity, aby
zmienić kąt kamery i obejrzeć ruch z innej perspektywy.
