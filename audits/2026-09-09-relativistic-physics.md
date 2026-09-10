**Audyt relatywistycznej fizyki klasycznej — 9 września 2026**

Badana rewizja: `7889077`. Zakres: kinematyka, pola retardowane, siły i precesja dipoli, reakcja promieniowania, bilanse, historia źródeł, opcjonalne ZPF oraz relacja integratora mechanicznego do estymatora sekularnego. Przegląd kodu połączono z istniejącą walidacją i niezależnymi sondami analitycznymi. Nie zmieniano równań ani parametrów produkcyjnych. Nie wykonywano pełnego skanu wszystkich par, energii i trybów statystycznych.

**Werdykt:** projekt zawiera poprawne składniki relatywistycznej elektrodynamiki, ale nie stanowi kompletnego, spójnego we wszystkich układach odniesienia modelu klasycznych cząstek z dipolami. Potwierdzono trzy defekty. Najpoważniejszy dotyczy transformacji siły dipolowej i występuje już przy umiarkowanych prędkościach. Wyników tych nie należy przenosić liczbowo na czasy kolapsu pozytonium bez osobnego pomiaru: podane procenty odnoszą się do izolowanych sond.

**1. [P1] Siła `covariantDipoleGradientForce` nie jest ogólnie kowariantna.**

Miejsca: `modules/electrodynamics.hpp:2925–2926,3174`; geometria dotychczasowego testu: `modules/maxwell_validation.hpp:1595–1604,2410–2429`.

Funkcja różniczkuje sprzężenie `U = mu_lab·B + p_lab·E`, po czym zwraca `grad(U)/gamma`. Równanie ruchu traktuje wynik jako `dp/dt` dla stałej masy. Samo wydzielenie części przestrzennej gradientu skalaru i podzielenie przez gamma nie zapewnia zgodności z tak zdefiniowaną czterosiłą: składowa czasowa fizycznej siły jest związana z pracą, a gradientu — z pochodną czasową pola.

Sonda: nieruchomy dipol źródłowy w początku układu, moment w osi z; dipol próbny w `x=2a0`, poruszający się wzdłuż x. Ładunek źródła wyzerowano wyłącznie w sondzie, aby odizolować sektor dipolowy. Historie są zadanymi prostoliniowymi liniami świata; oba opisy dotyczą tego samego zdarzenia obserwacji z uwzględnieniem względności jednoczesności. Boost jest również wzdłuż x. Dla siły i prędkości współliniowych transformacja relatywistyczna wymaga `F'x = Fx`.

| Prędkość dipola próbnego/c | Prędkość boostu/c | Zmierzony F'x/Fx | Wartość wymagana |
|---|---|---|---|
| 0 | 0,35 | 1,00000834 | 1 |
| 0,2 | 0,35 | 1,07527779 | 1 |
| 0,2 | 0,6 | 1,13637311 | 1 |
| 0,6 | 0,35 | 1,26583334 | 1 |
| 0,6 | 0,6 | 1,56251303 | 1 |

Wynik odpowiada błędnemu czynnikowi `1/(1−Vv/c²)`; kontrola ze spoczywającym celem daje tło około `8,3e−6`, związane z numeryką/regularyzacją. Efekt rzędu 7,5–56,3% jest od niego wielokrotnie większy. Poprawienie znaku sprzężenia, opisane wcześniej w kodzie, nie rozwiązało całego problemu.

Istniejąca bramka daje `4,50856e−6`, lecz korzysta z orbity o prędkościach w osi x, separacji w osi y i boostu w osi z. Nie zastępuje próby z niezerową pracą siły w kierunku boostu. To defekt równania, którego ta geometria nie wykrywa.

Zalecenie: wyprowadzić wspólną kowariantną dynamikę siły, spinu i energii wewnętrznej przy jawnej konwencji pędu. Nie dopasowywać współczynnika poprawki do jednej sondy. Dodać próby współliniowe oraz pola zmienne w czasie, badające wszystkie cztery składowe. Kowariantne rozszerzenia siły Stern–Gerlacha zawierają dodatkowe pochodne, a ich zgodność z dynamiką spinu wymaga osobnego sprawdzenia; patrz [Rafelski, Formanek i Steinmetz, 2018](https://doi.org/10.1140/epjc/s10052-017-5493-2).

**2. [P2] Pole dipolowe zmienia fizyczną prędkość poprawnych źródeł powyżej 0,999c.**

Miejsca: `modules/electrodynamics.hpp:1763,1985`. Obcięcie działa w konstrukcji pola bliskiego i dalekiego. Komentarz uzasadnia je ochroną przed niefizycznym ruchem fikcyjnych biegunów, lecz warunek obejmuje również poprawny jednostajny ruch rzeczywistego źródła. Położenie retardowane i moment nadal odpowiadają pierwotnej prędkości, natomiast pole jest liczone z innej prędkości.

Niezależne odniesienie: analityczne pole statycznego dipola przetransformowane do układu laboratoryjnego. Źródło porusza się wzdłuż y, moment jest w osi z, obserwacja w `(2a0,0,0)` w chwili zero. Wtedy `Bz = −gamma*mu0*m0/(4*pi*R³)`. Sonda ma stały moment i brak fizycznego przyspieszenia, więc ogranicznik nie ma tu uzasadnienia.

| beta źródła | Bz z kodu / Bz odniesienia |
|---|---|
| 0,3 | 1,0000000005 |
| 0,8 | 1,0000000004 |
| 0,999 | 0,999999982 |
| 0,9995 | 0,592592558 |
| 0,9999 | 0,060085498 |

Zaniżenie wynosi odpowiednio około 40,7% i 94,0% w dwóch ostatnich przypadkach. Jest to problem zakresu ultrarelatywistycznego, nie wykazany błąd powolnej orbity pozytonium.

Zalecenie: odróżnić nieważną kinematykę biegunów od poprawnej prędkości źródła; poprawić rekonstrukcję albo jawnie odrzucić stan spoza obsługiwanego zakresu. Wzorzec testowy musi pozostać niezależnym boostem pola statycznego. Dzisiejsze porównanie dwóch separacji tych samych fikcyjnych biegunów (`maxwell_validation.hpp:2852–2920`) mierzy zbieżność konstrukcji, ale wspólnego błędu fizycznego obu obliczeń nie wykryje. [Sautbekov, 2018](https://arxiv.org/html/1806.07089) stanowi dodatkowy punkt odniesienia dla pól i promieniowania ruchomych dipoli; sama zgodność współczynnika transformacji momentu nie certyfikuje całego solvera.

**3. [P2, opcjonalne `--zpf`] Pole ze zmiennym pasmem nie spełnia prawa Faradaya.**

Miejsce: `modules/zero_point_field.hpp:70–89`. Częstotliwość i amplituda zależą od aktualnej orbity: `k(t)=f*omega(t)/c`, `A(t)∝omega(t)²`, faza czasowa jest całką częstotliwości. Jednocześnie funkcja przyjmuje `B=n×E/c` jak dla monochromatycznej fali płaskiej. Te założenia nie tworzą rozwiązania Maxwella przy zmiennym omega.

Wystarczy jeden mod: `n=z`, `E` w osi x, `B` w osi y. W początku układu i w maksimum cosinusa `curl(E)=0`, natomiast `dB/dt=A'(t)/c` jest niezerowe. Zatem `dB/dt+curl(E) != 0`. Całkowanie fazy usuwa skoki, lecz nie usuwa tego członu ani pochodnych zmiennego wektora falowego.

Sonda używa `omega(t)=10^16*(1+10^14*t) s^−1`, `A(0)=1 V/m`. Zmniejszanie kroku `10^−18 → 10^−19 → 10^−20 s` daje residuum Faradaya podzielone przez analityczne `2*10^14/c`: `0,999925 → 0,99999925 → 0,9999999925`. Błąd zbiega do wartości niezerowej; nie wynika z różniczkowania numerycznego.

Zalecenie: traktować obecną konstrukcję jako zadane wymuszenie fenomenologiczne albo przebudować ją na pole z potencjałów spełniających Maxwella. Wynik tego wariantu nie jest rozstrzygającym testem rzeczywistej próżni SED. Pełne widmo proporcjonalne do omega³ i skończone pasmo zależne od trajektorii mają różne własności. Brak nie dotyczy przebiegów z wyłączonym ZPF; nie kwestionuje też poprawności pojedynczego modu o stałej częstotliwości.

**Kompletność i elementy, w których nie stwierdzono nowego defektu**

| Sektor | Ocena i granica wniosku |
|---|---|
| Kinematyka | Relatywistyczne p=gamma*m*v, energia kinetyczna, odwrotna konwersja p→v i kinematyka dwuciałowa są obecne. Nie znaleziono błędu wzorów w badanym zakresie. |
| Ładunki | Pełne człony prędkościowe i przyspieszeniowe pola Liénarda–Wiecherta oraz siła Lorentza. Historia i retardacja są przybliżane numerycznie. |
| Precesja BMT | Poprawna standardowa postać dla momentu spoczynkowego i dokładny obrót przy zamrożonym polu/prędkości. Nie dowodzi to poprawności całej sprzężonej dynamiki z siłą gradientową i reakcją promieniowania. |
| Reakcja ładunków | Zredukowana reakcja LL i projekcja ortogonalna czterosiły. Uniknięto ponownego dodania wzajemnej reakcji, gdy partner już wytwarza pole retardowane. LL pozostaje przybliżeniem, nie pełnym rozwiązaniem problemu samooddziaływania. |
| Promieniowanie M1/E2 | Są człony multipolowe i koherencja, lecz lokalna reakcja M1 jest długofalowa (`electrodynamics.hpp:1217–1274`). Nie stanowi pełnej relatywistycznej reakcji dowolnie przyspieszającego dipola. |
| Energia i pęd | Energia/pędy Darwina są przybliżone do rzędu v²/c². Nie są dokładnymi całkami pełnej retardowanej dynamiki. `boundField*` jest resztą domykającą bilans (`electrodynamics.hpp:3799`), a nie niezależnym pomiarem energii/pędu pola. To ograniczenie jest już wyraźnie opisane w README. |
| Historia i regulator | Historia jest skończona; wcześniejsze stany są ekstrapolowane (`retarded_charge_kinematics.hpp:157–168`). Retencja rzędu 4r/c nie wystarcza uniwersalnie przy dużych boostach. Regulator zależy od odległości laboratoryjnej i w rdzeniu nie jest pełnym kowariantnym modelem rozciągłego źródła. |
| Błąd całkowania | Kontroler porównuje położenia i prędkości (`crem_engine.hpp:161–179`), pomijając błąd orientacji spinu, fazy ZPF i strumieni. Mały błąd orbity nie certyfikuje dokładności tych obserwabli. |
| Kolaps i emisja | Estymator sekularny używa orbit Keplera i uśrednionej mocy Larmora (`crem_collapse.hpp:433`); domyślna emisja korzysta z kwantów (`crem_trajectory.hpp:69`). To model hybrydowy, a nie bezpośrednie całkowanie kompletnego klasycznego układu Maxwell–cząstki. |

Redukcja rzędu LL ma określony zakres stosowalności; sam fakt użycia relatywistycznej notacji go nie rozszerza. Patrz [Ekman, Heinzl i Ilderton, 2021](https://arxiv.org/abs/2105.01640). W tym audycie nie utożsamiano krótkoczasowego residuum pracy reakcji z łamaniem zachowania energii: trzeba uwzględnić odwracalny człon Schotta i dokładność całkowania.

Brak anihilacji, kwantowej statystyki spinu oraz wyprowadzenia czasów życia para/orto nie jest brakiem w klasycznej relatywistycznej teorii: to granica jej zastosowania do rzeczywistego pozytonium. Momenty magnetyczne, g i warunki przygotowania mogą być wejściami empirycznymi, ale wtedy nie są wynikami klasycznego modelu.

**Wykonana weryfikacja i odtworzenie**

`./positronium_validation --statistics-profile small`: **53/53 PASS**, 46,13 s. Sprawdzono aktualność istniejącego pliku wykonywalnego względem źródeł przez `make -q validation`. Nie przebudowywano głównego programu. Wynik oznacza przejście zadanych bramek; niezależne sondy powyżej pokazują ich luki.

Ważne wyniki: normalizacja Larmora w produkcyjnym próbkowaniu `1,0000456`, akumulacja `0,99997216`; dryf normy BMT `1,33e−15`; niezależny dwuperiodowy bilans energii: residuum względem energii wypromieniowanej `5,94%`, po zaostrzeniu tolerancji `0,981%`. Ten ostatni spadek świadczy o istotnym udziale błędu numerycznego. Sonda bilansowa inicjalizuje zerowe dipole, więc nie certyfikuje wykrytego wyżej sektora dipolowego.

Kod sond: [tools/audit_relativistic_physics.cpp](../tools/audit_relativistic_physics.cpp). Sondy korzystają bezpośrednio z produkcyjnych nagłówków, bez ROOT. Mają charakter pomiarowy i wypisują wyniki; ich kod wyjścia nie jest bramką poprawności fizycznej.

```sh
g++ -std=c++20 -O2 -I . tools/audit_relativistic_physics.cpp -o /tmp/positronium-physics-audit
/tmp/positronium-physics-audit
```

Archiwum: [wyniki sond](2026-09-09-physics-probes.txt), [pełny log walidacji](2026-09-09-validation.txt). Priorytet dalszych prac: poprawne równanie siły dipolowej i niezależne testy boostu; następnie obsługa szybkich źródeł oraz osobna korekta opcjonalnego ZPF.
