# Klasyczne pozytonium — symulacja elektrodynamiczna

Program w CERN ROOT przedstawia klasyczny model elektrodynamiczny układu
elektron–pozyton:

- przyciąganie cząstek jest liczone w kontrolowanym, niskoprędkościowym
  przybliżeniu Coulomba (początkowo `v/c ≈ 0,005`);
- ruch jest całkowany dla relatywistycznego pędu `p = γmv`;
- reakcja promieniowania wynika z nierelatywistycznej redukcji rzędu
  Landaua–Lifshitza `F_RR = τ dF_C/dt`, z analityczną pochodną siły Coulomba;
- moc promieniowania jest niezależnie liczona wzorem Liénarda i służy do
  kontroli wyniku, ale nie jest używana do korygowania trajektorii.

Każda cząstka ma też klasyczny moment dipolowy o losowym kierunku początkowym.
Losowane są również radialna i styczna prędkość względna. Kontrolowane zakresy
losowania pozwalają uzyskać bezpośrednie zderzenie, rozproszenie albo związany
stan para-/orto-podobny. Rzeczywista precesja jest na tej skali czasu bardzo
wolna, dlatego orientacja strzałek może zmieniać się tylko nieznacznie. Momenty
dipolowe nie są włączone do ruchu orbitalnego: ich klasyczne oddziaływanie jest
tutaj wiele rzędów wielkości słabsze od przyciągania Coulomba, a spin
pozytonium jest w istocie zjawiskiem kwantowym.

To nadal jest przybliżony model punktowych ładunków: dipole reprezentują
fenomenologicznie spin, a redukcja Landaua–Lifshitza nie zastępuje pełnego
problemu samopola. Stabilność i anihilację prawdziwego pozytonium
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

Nad tabelą widoczny jest wynik klasycznego drzewa decyzyjnego oraz cztery
decydujące warunki początkowe: energia względna `E_rel`, orbitalny moment pędu
`L_orb`, przewidywana odległość największego zbliżenia `r_min` i czas życia
`Lifetime`. Dla rozproszenia czas życia jest nieskończony.

Przyciski nad widokiem: `STOP` wstrzymuje animację (zmienia się na `START`), a
`EXIT` natychmiast zamyka program.

Bilans energii (łącznie z energią bliskiego pola Schotta) i zakres odległości
można sprawdzić bez uruchamiania okna:

```bash
./positronium --diagnose
```

Konkretny losowy przypadek można odtworzyć przez podanie ziarna, np.:

```bash
./positronium --diagnose --seed 8
```

Symulacja jest renderowana w 3D. Przeciągnij myszą w obszarze orbity, aby
zmienić kąt kamery i obejrzeć ruch z innej perspektywy.
