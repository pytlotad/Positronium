# Positronium — trwałe instrukcje projektu dla PRISM

## 1. Zastosowanie

Ten plik zawiera nadrzędne instrukcje korzystania z materiałów projektu
Positronium podczas przygotowywania publikacji naukowej w PRISM.

Kanoniczne repozytorium projektu:

https://github.com/pytlotad/Positronium

PRISM ma korzystać przede wszystkim z plików faktycznie zaimportowanych do
bieżącego projektu. Sam adres GitHub nie oznacza, że PRISM ma dostęp do
najnowszej zawartości repozytorium ani że projekt jest z nim synchronizowany.

## 2. Obowiązkowa kolejność pracy

Przed napisaniem, poprawieniem lub oceną publikacji PRISM ma:

1. przeczytać w całości ten plik;
2. przeczytać w całości `README.ap`, który jest głównym indeksem kodu, fizyki,
   ścieżek wykonania i konwencji nazewniczej wyników;
3. przeczytać odpowiednie części `README.md` opisujące równania, założenia,
   walidację i ograniczenia modelu;
4. sprawdzić implementację omawianego zagadnienia w `positronium.cpp` i w
   używanych przez niego plikach nagłówkowych;
5. sprawdzić `ScientificalReferences.txt` i kompletność każdej pozycji przed
   użyciem jej w bibliografii;
6. zinwentaryzować odpowiednie pliki z `distributions/` i interpretować ich
   nazwy wyłącznie według sekcji 8 pliku `README.ap`;
7. przed rozpoczęciem redagowania przedstawić krótki audyt wykorzystanych
   materiałów.

Audyt ma zawierać ścieżkę każdego wykorzystanego pliku, jego rolę, informację,
czy został odczytany, oraz listę brakujących lub niejednoznacznych danych.

## 3. Hierarchia źródeł wewnątrz projektu

W przypadku rozbieżności dotyczącej aktualnego działania programu obowiązuje
hierarchia opisana w `README.ap`: kod źródłowy, `README.ap`, `README.md`,
istniejące binaria, wygenerowane PDF-y. Binaria są artefaktami kompilacji, a
nie źródłem programu.

Dokumentacja i kod opisują model oraz implementację. Nie zastępują literatury
naukowej jako podstawy zewnętrznych twierdzeń fizycznych. Wykresy programu są
wynikami symulacji, a nie wynikami pomiarów.

## 4. Zasady interpretacji naukowej

PRISM ma zawsze:

- rozróżniać tryby `Visual simulation` i `Statistical analysis`;
- rozróżniać klasyczny model trajektorii od fenomenologicznego generatora
  rozpadów pozytonium;
- oddzielać założenia modelu, implementację numeryczną, wyniki symulacji,
  interpretację autora i wiedzę zaczerpniętą z literatury;
- nie przedstawiać modelu CREM jako kompletnego opisu QED;
- nie przedstawiać kanału `short-range` jako przekroju anihilacji QED;
- nie przedstawiać klasyfikacji `para` i `ortho` w trybie Visual jako
  kwantowych stanów singletowego i trypletowego;
- nie wyprowadzać niepodanych parametrów uruchomienia z samej nazwy PDF;
- oznaczać brak informacji albo niepewność zamiast wymyślać dane, źródła,
  parametry lub wnioski;
- wskazywać ścieżkę pliku projektu stanowiącego podstawę opisu implementacji;
- sprawdzać, czy podpis każdego użytego wykresu opisuje tryb, eksperyment,
  przedstawianą wielkość oraz kontekst utworzenia pliku.

Każdy PDF z `distributions/` należy traktować jako wynik ostatniego udanego
uruchomienia zapisującego daną nazwę. Pliki o różnych nazwach nie muszą
pochodzić z tej samej sesji. Liczbę zdarzeń, ziarno, energię wiązki i inne
parametry wolno podać tylko wtedy, gdy wynikają z adnotacji wykresu, zapisanego
polecenia uruchomienia albo innego jednoznacznego materiału projektu.

## 5. Atrybucja i dostępność

Program może być używany przez każdego. Publikacja wyników opartych na tym
programie lub na modelu CREM — także po zmianie kodu albo modelu — powinna
wskazywać:

- repozytorium: https://github.com/pytlotad/Positronium
- autora modelu CREM: Tadeusz Sławomir Pytlos
- kontakt: tadeusz.slawomir.pytlos@gmail.com

PRISM ma zachować tę informację w publikacji, materiałach uzupełniających i
opisie dostępności kodu. Nie wolno rozszerzać tego zalecenia o warunki prawne,
których nie podano w repozytorium.

Zalecany tekst sekcji dostępności kodu w publikacji:

> The source code, CREM model documentation, simulation methodology, and
> generated graphical materials used in this work are available at
> https://github.com/pytlotad/Positronium. The CREM model was developed by
> Tadeusz Sławomir Pytlos (tadeusz.slawomir.pytlos@gmail.com).

## 6. Kontrola wersji materiałów

Jeżeli w projekcie publikacji zapisano identyfikator commitu, PRISM ma go
podawać jako wersję analizowanych materiałów. Jeżeli identyfikatora brakuje,
PRISM ma napisać, że analizuje zaimportowaną migawkę o nieustalonej wersji i
poprosić o identyfikator commitu przed sformułowaniem ostatecznej informacji o
reprodukowalności.

PRISM nie może twierdzić, że sprawdził aktualny stan GitHuba, jeżeli korzystał
wyłącznie z plików znajdujących się w projekcie.

## 7. Odpowiedź inicjalizacyjna PRISM

Po otrzymaniu polecenia rozpoczęcia pracy PRISM ma najpierw odpowiedzieć:

1. czy odczytał `PRISM_INSTRUCTIONS.md` i `README.ap` w całości;
2. jakie pliki źródłowe, dokumentacyjne, bibliograficzne i graficzne widzi;
3. jaką wersję lub commit materiałów potrafi potwierdzić;
4. jakich materiałów brakuje;
5. dopiero potem — jaki plan publikacji proponuje.

Jeżeli PRISM nie może odczytać któregoś obowiązkowego pliku, ma zatrzymać
redagowanie i jednoznacznie wskazać brakujący plik.
