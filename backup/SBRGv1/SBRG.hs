-- TODO: below didn't work :(
-- TODO go back to save restore necessary trash before increasing max_terms
-- This is done by preferring terms that will be closer to having k=0,3 for the used_is near the end of the RG. This should scale as O(n) if we apply (and later unapply) c4 to the RG data

-- ghci -fbreak-on-error
-- :set -fbreak-on-error

{-# LANGUAGE TupleSections #-}
-- :set -XTupleSections

{-# OPTIONS -Wall -fno-warn-unused-binds -fno-warn-unused-imports -O2 -optc-O2 -optc-march=native -optc-mfpmath=sse #-}
-- -rtsopts -prof -fprof-auto        -ddump-simpl -threaded
-- +RTS -N4 -xc -sstderr -p

import Control.Monad
import qualified Control.Monad.ST as ST
import Data.Bits (testBit, xor)
import Data.Either
import Data.Function
import Data.List
import Data.Maybe
import Data.Tuple
import Debug.Trace
import System.Environment
import System.CPUTime

--import qualified Control.Parallel.Strategies as Parellel

--import Control.Monad.Trans.State.Strict (State)
import qualified Control.Monad.Trans.State.Strict as State

import Data.IntSet (IntSet)
import qualified Data.IntSet as IntSet

import Data.Set (Set)
import qualified Data.Set as Set hiding (fromList)

import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IntMap hiding (fromList, insert, delete, adjust, adjustWithKey, update, updateWithKey) -- use alter instead

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map hiding (fromList, insert, delete, adjust, adjustWithKey, update, updateWithKey) -- use alter instead

-- import Data.Sequence (Seq)
-- import qualified Data.Sequence as Seq

import qualified Data.Vector as Vector
import qualified System.Random.MWC as Random
import qualified System.Random.MWC.Distributions as Random

toF :: Double -> F

-- type F = LogFloat
-- toF = fromDouble'LF

type F = Double
toF = id


-- debug

debug :: Bool
debug = True

check :: Bool -> String -> a -> a
check b s x = not debug || b ? x $ error s

-- generic

if' :: Bool -> a -> a -> a
if' True  x _ = x
if' False _ y = y

infixr 1 ?
(?) :: Bool -> a -> a -> a
(?) = if'

boole :: Num a => Bool -> a
boole False = 0
boole True  = 1

(.*) :: (c -> c') -> (a -> b -> c) -> a -> b -> c'
(.*) = (.) . (.)

fromLeft :: Either a b -> a
fromLeft (Left x) = x
fromLeft _        = error "fromLeft"

fromRight :: Either a b -> b
fromRight (Right y) = y
fromRight _         = error "fromRight"

fst3 :: (a,b,c) -> a
fst3 (x,_,_) = x

snd3 :: (a,b,c) -> b
snd3 (_,y,_) = y

thd3 :: (a,b,c) -> c
thd3 (_,_,z) = z

mapPair :: (a -> c, b -> d) -> (a,b) -> (c,d)
mapPair (f,g) (x,y) = (f x, g y)

mapBoth :: (a -> b) -> (a,a) -> (b,b)
mapBoth f = mapPair (f,f)

mapFst :: (a -> c) -> (a,b) -> (c,b)
mapFst f = mapPair (f, id)

mapSnd :: (b -> c) -> (a,b) -> (a,c)
mapSnd f = mapPair (id, f)

-- getArgs3 :: (Read a, Read b, Read c) => IO (a,b,c)
-- getArgs3 = do
--   [x1,x2,x3] <- getArgs
--   return (read x1, read x2, read x3)

-- getArgs4 :: (Read a, Read b, Read c, Read d) => IO (a,b,c,d)
-- getArgs4 = do
--   [x1,x2,x3,x4] <- getArgs
--   return (read x1, read x2, read x3, read x4)

getArgs5 :: (Read t1, Read t2, Read t3, Read t4, Read t5) => IO (t1,t2,t3,t4,t5)
getArgs5 = do
  [x1,x2,x3,x4,x5] <- getArgs
  return (read x1, read x2, read x3, read x4, read x5)

infixl 1 `applyIf`
applyIf :: (a -> a) -> (a -> Bool) -> a -> a
applyIf f q x = if' (q x) f id $ x

nest :: Int -> (a -> a) -> a -> a
nest 0 _ = id
nest n f = nest (n-1) f . f

justIf :: (a -> Bool) -> a -> Maybe a
justIf q x = q x ? Just x $ Nothing

meanError :: Floating f => [f] -> (f,f)
meanError xs = (mean, sqrt $ (x2/x0 - mean*mean)/(x0-1))
  where [x0,x1,x2] = map (\n -> sum $ map (^n) xs) [0,1,2::Int]
        mean = x1/x0

epsilon, epsilon_2 :: Fractional f => f
epsilon   = 2 ^^ (-52::Int)
epsilon_2 = 2 ^^ (-26::Int)

-- Nothing if x+y is tiny compared to x and y
infixl 6 +?
(+?) :: F -> F -> Maybe F
x +? y = x_y*x_y/(x*x + y*y) < epsilon ? Nothing $ Just x_y -- TODO error "+? actually filtered something!"
  where x_y = x + y

insertAdd'Map :: Ord k => k -> F -> Map k F -> Map k F
insertAdd'Map kx x = Map.alter (maybe (Just x) (+?x)) kx -- Map.insertWith (+) kx x

unionAdd'Map :: Ord k => Map k F -> Map k F -> Map k F
unionAdd'Map m1 m2 = Map.mergeWithKey (\_ x y -> x +? y) id id m1 m2

unionsAdd'Map :: Ord k => [Map k F] -> Map k F
unionsAdd'Map = foldl' unionAdd'Map Map.empty

-- LogFloat
-- https://hackage.haskell.org/package/logfloat-0.13.3.1/docs/Data-Number-LogFloat.html#v:product

data LogFloat = LogFloat Bool Double
  deriving (Eq)

toDouble'LF :: LogFloat -> Double
toDouble'LF (LogFloat s y) = s ? exp y $ -exp y

fromDouble'LF :: Double -> LogFloat
fromDouble'LF x = LogFloat (x>=0) (log $ abs x)

instance Show LogFloat where
  showsPrec n z@(LogFloat s y) str
    | abs e > 300 && not (isInfinite y || isNaN y)
                = if' s "" "-" ++ shows (exp $ ln10*y') ('e' : shows (e::Integer) str)
    | otherwise = showsPrec n (toDouble'LF z) str
    where (e,y') = properFraction $ y / ln10
          ln10   = log 10

instance Read LogFloat where
  readsPrec n = map (mapFst fromDouble'LF) . readsPrec n

instance Ord LogFloat where
  compare (LogFloat s1 ln1) (LogFloat s2 ln2)
    | not s1 &&     s2 = LT
    |     s1 && not s2 = GT
    |     s1 &&     s2 = compare ln1 ln2
    | otherwise        = compare ln2 ln1

foreign import ccall unsafe "math.h log1p"
    log1p :: Double -> Double

instance Num LogFloat where
  z1@(LogFloat s1 y1) + z2@(LogFloat s2 y2)
    | isInfinite y1 && z1==z2 = z1
    | y1 >= y2                = LogFloat s1 $ y1 + log1p (if' (s1==s2) id negate . exp $ y2 - y1)
    | otherwise               = z2 + z1
  (LogFloat s1 y1) * (LogFloat s2 y2) = LogFloat (not $ xor s1 s2) (y1 + y2)
  abs (LogFloat _ y) = LogFloat True y
  signum z@(LogFloat s y)    = isInfinite y && y < 0 ? z $ LogFloat s 0
  fromInteger = fromDouble'LF . fromInteger
  negate (LogFloat s y) = LogFloat (not s) y

instance Fractional LogFloat where
  fromRational = fromDouble'LF . fromRational
  recip (LogFloat s y) = LogFloat s (-y)

instance Floating LogFloat where
  sqrt (LogFloat True y) = LogFloat True (y/2)
  sqrt _                 = error "LogFloat sqrt"
  pi    = undefined
  exp   = undefined
  log   = undefined
  sin   = undefined
  cos   = undefined
  asin  = undefined
  acos  = undefined
  atan  = undefined
  sinh  = undefined
  cosh  = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined

-- rankZ2

rankZ2 :: [[Bool]] -> Int
rankZ2 = go 0 . map (foldl' (\x b -> 2*x + boole b) 0)
  where
    go :: Int -> [Integer] -> Int
    go n []         = n
    go n (0:rows)   = go n rows
    go n (row:rows) = go (n+1) $ map (xor row `applyIf` flip testBit j) rows
      where
        j = head $ filter (testBit row) [0..]

-- Sigma

type Sigma = IntMap Int
type SigmaTerm = (Sigma, F)

check'Sigma :: Sigma -> Bool
check'Sigma g = all (\k -> 1<=k && k<=3) $ IntMap.elems g

toList'GT :: SigmaTerm -> ([(Int,Int)],F)
toList'GT = mapFst IntMap.toList

fromList'GT :: ([(Int,Int)],F) -> SigmaTerm
fromList'GT = mapFst $ IntMap.fromListWith $ error "fromList'GT"

toLists'mapGT :: Map Sigma F -> [([(Int,Int)],F)]
toLists'mapGT = map toList'GT . Map.toList

fromList'mapGT :: [SigmaTerm] -> Map Sigma F
fromList'mapGT = Map.fromListWith (+)

fromLists'mapGT :: [([(Int,Int)],F)] -> Map Sigma F
fromLists'mapGT = fromList'mapGT . map fromList'GT

show'mapGT :: Map Sigma F -> String
show'mapGT gTs = "fromLists'mapGT " ++ show (toLists'mapGT gTs)

scaleST :: F -> SigmaTerm -> SigmaTerm
scaleST x (g,h) = (g,x*h)

acommQ :: Sigma -> Sigma -> Bool
acommQ g1 g2 = IntMap.foldl' xor False $ IntMap.intersectionWith (/=) g1 g2

-- the first returned Int is the phase in units of pi/2
multSigma1 :: Int -> Int -> (Int,Int)
multSigma1 1 1 = ( 0,0)
multSigma1 2 2 = ( 0,0)
multSigma1 3 3 = ( 0,0)
multSigma1 1 2 = ( 1,3)
multSigma1 2 1 = (-1,3)
multSigma1 2 3 = ( 1,1)
multSigma1 3 2 = (-1,1)
multSigma1 3 1 = ( 1,2)
multSigma1 1 3 = (-1,2)
multSigma1 k 0 = ( 0,k)
multSigma1 0 k = ( 0,k)
multSigma1 _ _ = error "multSigma1"

multSigma :: Sigma -> Sigma -> (Int,Sigma)
multSigma g1 g2 = (IntMap.foldl' (+) 0 $ IntMap.intersectionWith (fst .* multSigma1) g1 g2,
                   IntMap.filter (/=0) $ IntMap.unionWith (snd .* multSigma1) g1 g2)

-- i [a,b]/2
icomm :: Sigma -> Sigma -> Maybe SigmaTerm
icomm g1 g2 = case mod n 4 of
                   0 -> Nothing
                   2 -> Nothing
                   1 -> Just (g,-1)
                   3 -> Just (g, 1)
                   _ -> error "icomm"
  where (n,g) = multSigma g1 g2

-- {a,b}/2
acomm :: Sigma -> Sigma -> Maybe SigmaTerm
acomm g1 g2 = case mod n 4 of
                   0 -> Just (g, 1)
                   2 -> Just (g,-1)
                   1 -> Nothing
                   3 -> Nothing
                   _ -> error "acomm"
  where (n,g) = multSigma g1 g2

-- c4 g4 g = R^dg g R where R = (1 + i g4)/sqrt(2)
c4 :: Sigma -> Sigma -> Maybe SigmaTerm
c4 g4 g = icomm g g4

-- TODO very slow
c4s :: [Sigma] -> SigmaTerm -> SigmaTerm
c4s g4s gT = foldl (\gT'@(g,c) g4 -> maybe gT' (scaleST c) $ c4 g4 g) gT g4s

local'Sigma :: Int -> Sigma -> Bool
local'Sigma n = \g -> IntMap.size g < sqrt_n
  where sqrt_n = round $ sqrt $ (fromIntegral n :: Double)

-- Ham

-- TODO make cg a maybe?
data Ham = Ham {
  gc'Ham   :: Map Sigma F,        -- Hamiltonian:  sigma matrix  -> coefficient
  cg'Ham   :: Map F (Set Sigma),  -- ordering:     |coefficient| -> sigma matrix
  igs'Ham  :: IntMap (Set Sigma), -- local sigmas: site index    -> sigma matrices acting on that site
  nlgs'Ham :: Set Sigma,          -- non-local sigmas
  ls'Ham   :: [Int] }             -- system lengths
  deriving (Eq)

n'Ham :: Ham -> Int
n'Ham = product . ls'Ham

instance Show Ham where
  showsPrec _ ham s = "fromLists'Ham " ++ show (ls'Ham ham) ++ " " ++ showsPrec 10 (toLists'Ham ham) s

check'Ham :: Ham -> Bool
check'Ham ham = ham == (fromLists'Ham (ls'Ham ham) $ toLists'Ham ham) && all check'Sigma (Map.keys $ gc'Ham ham)

zero'Ham :: [Int] -> Ham
zero'Ham ls = Ham Map.empty Map.empty IntMap.empty Set.empty ls

null'Ham :: Ham -> Bool
null'Ham = Map.null . gc'Ham

toLists'Ham :: Ham -> [([(Int,Int)],F)]
toLists'Ham = toLists'mapGT . gc'Ham

fromList'Ham :: [Int] -> [SigmaTerm] -> Ham
fromList'Ham ls = foldl' insert'Ham (zero'Ham ls) . Map.toList . fromList'mapGT

fromLists'Ham :: [Int] -> [([(Int,Int)],F)] -> Ham
fromLists'Ham ls = fromList'Ham ls . map fromList'GT

insert'Ham :: Ham -> SigmaTerm -> Ham
insert'Ham ham@(Ham gc cg igs nlgs ls) (g,h)
  | isNothing h_ =
    Ham (Map.insertWith (error "insert'Ham") g h gc)
        (Map.insertWith Set.union (abs h) (Set.singleton g) cg)
        (if' localQ igs' igs  )
        (if' localQ nlgs nlgs')
        ls
  | Just h' <- fromJust h_ +? h =
    Ham (Map.alter (const (Just h') . fromJust) g gc)
        (Map.insertWith Set.union (abs h') (Set.singleton g)
        $Map.alter (justIf (not . Set.null) . Set.delete g . fromJust) (abs $ fromJust h_) cg)
        igs nlgs ls
  | otherwise = fst $ deleteLookup'Ham ham g
  where n        = n'Ham ham
        h_       = Map.lookup g gc
        localQ   = local'Sigma n g
        igs'     = IntMap.unionWith Set.union igs $ IntMap.fromSet (const $ Set.singleton g) $ IntMap.keysSet g
        nlgs'    = Set.insert g nlgs

deleteLookup'Ham :: Ham -> Sigma -> (Ham, F)
deleteLookup'Ham ham@(Ham gc cg igs nlgs ls) g =
  (Ham gc'
       (alter' (justIf (not . Set.null) . Set.delete g . fromJust) (abs h) cg)
       (if' localQ igs' igs  )
       (if' localQ nlgs nlgs')
       ls, h)
  where n             = n'Ham ham
        (Just h, gc') = Map.updateLookupWithKey (\_ _ -> Nothing) g gc
        igs'          = foldl' (flip $ IntMap.alter $ justIf (not . Set.null) . Set.delete g . fromJust) igs $ IntMap.keys g
        nlgs'         = Set.delete g nlgs
        alter' a b c  = Map.alter a b c
        localQ        = ((==) `on` Set.size) nlgs nlgs'
      --localQ        = local'Sigma n g

-- union'Ham :: Ham -> Map Sigma F -> Ham
-- union'Ham ham gs0 = foldl' insert'Ham ham $ Map.toList gs0
-- 
-- difference'Ham :: Ham -> Set Sigma -> (Ham, Map Sigma F)
-- difference'Ham (Ham gc cg igs nlgs ls n) gs0 = (Ham gc' cg' igs' nlgs' ls n, gc0')
--   where gc0'  = Map.intersectionWith (curry fst) gc $ Map.fromSet (const ()) gs0
--         gc'   = Map.difference gc gc0'
--         cg'   = Map.differenceWith just_Set_difference cg $ Map.fromListWith Set.union $ map (\(g,c) -> (abs c,Set.singleton g)) $ Map.toList gc0'
--         igs'  = IntMap.differenceWith just_Set_difference igs $ IntMap.unionsWith Set.union
--               $ map (\g -> IntMap.fromSet (const $ Set.singleton g) $ IntMap.keysSet g) $ Set.toList $ Set.difference gs0 nlgs
--         nlgs' = Set.difference nlgs gs0
--         just_Set_difference x y = Just $ Set.difference x y

nearbySigmas :: Ham -> Sigma -> Set Sigma
nearbySigmas ham g = Set.union localSigmas nonlocalSigmas
  where localSigmas    = Set.unions $ IntMap.elems $ IntMap.intersectionWith (curry fst) (igs'Ham ham) g
        nonlocalSigmas = Set.filter (not . IntSet.null . (IntSet.intersection `on` IntMap.keysSet) g) $ nlgs'Ham ham

-- c4'Ham :: F -> Ham -> Sigma -> Ham
-- c4'Ham dir ham g4 = union'Ham ham' $ fromList'mapGT $ Map.elems $ Map.intersectionWith (\(g,c) c' -> (g,dir*c*c')) ggT gc0'
--   where -- map from old Sigma to new Sigma with a sign
--         ggT :: Map Sigma SigmaTerm
--         ggT = Map.mapMaybe id $ Map.fromSet (c4 g4) $ nearbySigmas ham g4
--         -- ham w/o old Sigmas and old Sigmas with coefficients
--         (ham', gc0') = difference'Ham ham $ Map.keysSet ggT

-- TODO 50% of CPU time is used in this function for MBL phases (TODO this may no longer be true)
-- perhaps using Map.union and Map.difference instead of insert and delete and using a monad for igs will help
c4'Ham :: F -> Sigma -> Ham -> Ham
c4'Ham dir g4 ham = foldl' insert'Ham ham' gTs'
  where (ham',gTs') = foldl' delSigma (ham,[]) $ Set.elems $ nearbySigmas ham g4
        delSigma (ham0,gTs) g = case c4 g4 g of
                                    Nothing      -> (ham0,gTs)
                                    Just (g',c') -> mapSnd (\c -> (g',dir*c*c'):gTs) $ deleteLookup'Ham ham0 g

-- TODO one could transpose (ie apply all c4 to one Sigma at a time instead of one c4 to all Sigma) this function by calculating the locality of the g4s
c4s'Ham :: F -> [Sigma] -> Ham -> Ham
c4s'Ham dir g4s ham = foldl' (flip $ c4'Ham dir) ham g4s

--

data RG = RG {
  ham'RG       :: Ham,
  diag'RG      :: Map Sigma F,
  unusedIs'RG  :: IntSet,
  g4s'RG       :: [Sigma],
  _unusedIs'RG :: IntSet,           -- future unusedIs
  _g4s'RG      :: [Sigma],          -- future reversed g4s
  trash'RG     :: [[SigmaTerm]],
  stab0'RG     :: [SigmaTerm],      -- stabilizers in new basis
  stab'RG      :: Ham,              -- stabilizers in old basis
  max_terms'RG :: Maybe Int}        -- used to lower insert'Ham.igs' time

ls'RG :: RG -> [Int]
ls'RG = ls'Ham . ham'RG

n'RG :: RG -> Int
n'RG = n'Ham . ham'RG

instance Show RG where
  showsPrec _ (RG ham diag unusedIs g4s _ _ trash _ stab _) s =
    "RG {ham'RG      =     "                 ++ show ham                     ++ ",\n" ++
    "    diag'RG     =     "                 ++ show'mapGT diag              ++ ",\n" ++
    "    unusedIs'RG =     IntSet."          ++ show unusedIs                ++ ",\n" ++
    "    g4s'RG      = map IntMap.fromList " ++ show (map IntMap.toList g4s) ++ ",\n" ++
    "    trash'RG    =                     " ++ show trash                   ++ ",\n" ++
    "    stab'RG     =     "                 ++ show stab                    ++ "}"   ++ s

init'RG :: Ham -> RG
init'RG ham = rg
  where rg = RG ham Map.empty (IntSet.fromList [0..(n'Ham ham - 1)]) [] IntSet.empty [] [] [] (stabilizers rg) (Just 64)
-- (\n -> Just $ 2*n) -- min 16 $ 

-- Everett's modified Ising model at B:(J=1,K=1,h=0) needs nonEpsilon=True for large Gamma due to exponentially small stabilizers
nonEpsilon :: F -> Bool
nonEpsilon _ = True
--nonEpsilon x = epsilon_2 < abs x -- epsilon_2

-- g ~ sigma matrix, _G ~ Sigma, i ~ site index, h ~ energy coefficient
rgStep :: RG -> RG
rgStep rg@(RG ham1 diag unusedIs g4s _unusedIs _g4s trash stab0 _ max_terms)
  | IntSet.null unusedIs = rg
  | null'Ham ham1        = let (i,unusedIs_)   = IntSet.deleteFindMin unusedIs
                               rg0 = rg {unusedIs'RG = unusedIs_, stab0'RG = (IntMap.singleton i 3,0):stab0, stab'RG = stabilizers rg0}
                           in  rg0
--   | null'Ham ham1 = rg
  | otherwise     = rg'
  where
    _G1            :: [(SigmaTerm,Sigma)]
    _G2            :: [(SigmaTerm,SigmaTerm)]
    _G3 ,_G4       :: Map Sigma F
    _G5 ,_G6, _G6' :: [SigmaTerm]
  --max_terms_     = Nothing -- Just 32 -- used to lower insertAdd'Map time. needed for eg Ising ls=[4,4] couplings=[0.1,0.1,1] Gamma=1 seed=0
  --max_terms0     = max_terms $ length _G2 -- used to lower insert'Ham.igs' time
    g3             = Set.findMin $ snd $ Map.findMax $ cg'Ham ham1
    i3             = fromMaybe (IntSet.findMin unusedIs) $ IntSet.lookupGE (fst $ IntMap.findMin g3) unusedIs
    g3'            = IntMap.singleton i3 3
    unusedIs'      = IntSet.delete i3 unusedIs
    g4s'           = find_G4s g3 g3'
                  -- apply c4 transformation, remove g3' from ham and extract its coefficient:
    (ham2,h3')     = flip deleteLookup'Ham g3' $ c4s'Ham 1 g4s' ham1
                  -- find sigma matrices that overlap with g3', split into anticommuting commutators (_G1) and diagonal matrices (diag')
    (diag0',_G1)   = partitionEithers $ mapMaybe (\g -> maybe (isDiag g ? Just (Left g) $ Nothing) (Just . Right . (,g)) (icomm g3' g))
                                      $ Set.toList $ nearbySigmas ham2 g3'
                  -- remove diagonal matrices from ham and extract coefficients for diag1'
    (ham3,diag1')  = foldl' (\(ham,gTs) g -> mapSnd (\h -> ((g,h):gTs)) $ deleteLookup'Ham ham g) (ham2,[]) diag0'
                  -- remove anticommuting matrices (_G2 ~ Delta) from ham and extract/calculate coefficients for _G2
    (ham4,_G2)     = foldl' strip_G (ham3,[]) _G1
                  -- keep only max_terms_ largest terms and distribute _G2
    _G3            = fst $ mult_G Map.empty _G2
                     -- $ maybe id (\_max_terms -> take _max_terms . sortOn (\((_,hL),(_,hR)) -> -abs (hL*hR))) max_terms_ $ _G2
                  -- do Everett's growth rate cut
  --_G3'           = Map.fromList $ take (2*length _G2) $ sortOn (negate . abs . snd) $ Map.toList _G3
                  -- extract diagonal terms and include g3' with the diagonal terms
    (diag'',_G4)   = mapFst (insertAdd'Map g3' h3') $ Map.partitionWithKey (\g _ -> isDiag g) _G3
    (_G5,trash')   = partition (nonEpsilon . snd) $ Map.toList _G4
                  -- filter tiny terms, keep only max_terms terms, and add non-diagonal terms back into Hamiltonian
    (_G6 ,trash_)  = maybe (,[]) (\_max_terms -> splitAt _max_terms . sortOn (negate . abs . snd)) max_terms _G5
    (_G6',trash'') = IntSet.null _unusedIs || isNothing max_terms ? ([],trash_)
                   $ splitAt (fromJust max_terms) $ sortOn (mapPair (treasurePenalty, negate . abs)) trash_
      where treasurePenalty g = IntMap.size $ IntMap.filter (/=3) $ (fst $ c4s _g4s (g,0)) IntMap.\\ IntMap.fromSet (const ()) _unusedIs
  --(_G6',trash'') = splitAt (fromMaybe 0 max_terms) $ sortOn (\(g,c) -> (-size g,abs c)) $ trash_
  --  where size g = IntMap.size $ IntMap.intersection g $ IntMap.fromSet (const ()) unusedIs
    ham5           = foldl' (\ham (g,h) -> insert'Ham ham (g,h)) ham4 $ _G6 ++ _G6'
    rg'            = RG ham5 (unionsAdd'Map [fromList'mapGT diag1',diag'',diag]) unusedIs' (reverse g4s' ++ g4s)
                       _unusedIs (null _g4s ? _g4s $ tail _g4s) ((trash''++trash'):trash) ((g3',h3'):stab0) (stabilizers rg') max_terms
    
    isDiag :: Sigma -> Bool
    isDiag g = not $ any (flip IntSet.member unusedIs') $ IntMap.keys g
    
    find_G4s :: Sigma -> Sigma -> [Sigma]
    find_G4s g0 g1 = maybe (g0 == g1 ? [] $ find_G4s g_ g1 ++ find_G4s g0 g_) (replicate 1 . fst) $ icomm g0 g1
      where
        (i0s,i1s) = mapBoth (IntSet.intersection unusedIs . IntMap.keysSet) (g0,g1)
        i01s      = IntSet.intersection i0s i1s
        is        = map IntSet.findMin $ IntSet.null i01s ? [i0s,i1s] $ [i01s]
        g_        = IntMap.fromListWith (error "find_G4s") $ [(i, head $ [1,2,3] \\ map (IntMap.findWithDefault 0 i) [g0,g1]) | i <- is] -- |
    
    strip_G :: (Ham, [(SigmaTerm,SigmaTerm)]) -> (SigmaTerm,Sigma) -> (Ham, [(SigmaTerm,SigmaTerm)])
    strip_G (ham_, _G_) (gT,g) = let x = -1/h3' in mapSnd (\h -> (scaleST (x*h) gT,(g,h)):_G_) $ deleteLookup'Ham ham_ g
                                        -- the extra minus is because we apply icomm twice
    
    -- TOOD maybe a union approach would work here?
    mult_G :: Map Sigma F -> [(SigmaTerm,SigmaTerm)] -> (Map Sigma F,[(SigmaTerm,SigmaTerm)])
    mult_G _G_' ((gTL@(gL,hL),(gR,hR)):_G_tail) = mult_G (foldl' mult _G_' $ (gTL,(gR,0.5*hR)):_G_tail) _G_tail
      where
        mult :: Map Sigma F -> (SigmaTerm,SigmaTerm) -> Map Sigma F
        mult _G0' (_,(gR0,hR0)) = maybe _G0' (\(g',h') -> insertAdd'Map g' (h'*hL*hR0) _G0') $ icomm gL gR0
    mult_G _G_' [] = (_G_',[])

stabilizers :: RG -> Ham
stabilizers rg = c4s'Ham (-1) (g4s'RG rg) $ fromList'Ham (ls'RG rg) $ stab0'RG rg
  where n = n'RG rg

-- stabilizers :: RG -> Ham
-- stabilizers rg = c4s'Ham (-1) (g4s'RG rg) $ fromList'Ham (ls'RG rg) $ [(g, Map.findWithDefault 0 g $ diag'RG rg) | i <- [0..n-1], let g = IntMap.fromList [(i,3)]] -- |
--   where n = n'RG rg

runRG :: RG -> RG
runRG = until (IntSet.null . unusedIs'RG) $ rgStep

runRG_safe :: RG -> RG
runRG_safe rg0 = fromRight $ runRG' False rg0{max_terms'RG = Just 16}
  where runRG' :: Bool -> RG -> Either ([Sigma], IntSet, Maybe Int) RG
        runRG' had_trash rg@(RG{unusedIs'RG = unusedIs})
          | rewind {-|-} = let max_terms' = if' (IntSet.null $ _unusedIs'RG rg) id (fmap (2*)) $ max_terms'RG rg
                           in  trace ("\n\n"++show unusedIs++"\n\n") $ Left (g4s'RG rg, unusedIs, max_terms')
          | null_unsedIs = Right rg -- |
          | is_ok        = _runRG'
          | otherwise    = let new        = fromLeft _runRG'
                               _unusedIs' = snd3 new
                               _g4s'      = reverse $ take (IntSet.size unusedIs - IntSet.size _unusedIs') $ fst3 new
                           in  runRG' had_trash rg{_g4s'RG = _g4s', _unusedIs'RG = _unusedIs', max_terms'RG = thd3 new}
          where null_unsedIs = IntSet.null unusedIs
                rewind       = null'Ham (ham'RG rg) && not null_unsedIs && had_trash
                rg'          = rgStep rg
                new_trash    = not $ null $ head $ trash'RG $ rg'
                _runRG'      = runRG' (had_trash || new_trash) rg'
                is_ok        = had_trash || not new_trash || isRight _runRG'

recoverHam :: RG -> Ham
recoverHam rg = c4s'Ham (-1) (g4s'RG rg) $ fromList'Ham (ls'RG rg) $ Map.toList $ diag'RG rg

-- return: 1 == RMS of <g>^2 over all eigenstates
-- Ham: the stabilizer Ham
rmsQ :: Sigma -> Ham -> Bool
rmsQ g0 stab = not $ any (acommQ g0) $ nearbySigmas stab g0

anderson_corr :: Int -> Sigma -> Ham -> (Double,Double)
anderson_corr x0 g0 stab | isJust gg = meanError [boole $ rmsQ (g_ i) stab | i <- [0..n-1]] -- |
                         | otherwise = (0,0)
  where
    n    = n'Ham stab
    ls   = ls'Ham stab
    lY   = product $ tail ls
    gg :: Maybe SigmaTerm -- g(0) g(x0)
    gg  = acomm g0 $ IntMap.mapKeys (flip mod n . (+ x0*lY)) g0
    g_ :: Int -> Sigma -- g(i) g(i+x0)
    g_ i = IntMap.mapKeys (i_xs ls . zipWith (+) (xs_i ls i) . xs_i ls) $ fst $ fromJust gg

ee_simple :: IntSet -> Ham -> Double
ee_simple = error "ee_simple"

-- entanglement entropy: arguments: stabilizer Ham and list of region sizes
ee1d :: [Int] -> Ham -> [(Int,Double,Double)]
ee1d l0s stab = [uncurry (l0,,) $ meanError [regionEE (l0*lY) i | i <- [0,lY..n-lY]] | l0 <- l0s]
  where
    n  = n'Ham stab
    lY = product $ tail $ ls'Ham stab
    -- entanglement entropy of the region [i..i+l-1]
    regionEE :: Int -> Int -> Double
    regionEE l i = (0.5*) $ fromIntegral $ rankZ2 $ acommMat regionStabs
      where regionStabs = [IntMap.filterWithKey (\j _ -> inRegionQ j) g
                          | g <- (union `on` flip (IntMap.findWithDefault []) cutStabs') i (mod'$i+l) ++ nonlocalStabs, -- |
                            cutByRegion l i g]
            -- is j in the region [i..i+l-1] ?
            inRegionQ :: Int -> Bool
            inRegionQ j = mod' (j-i) < l
    mod' i = mod i n
    acommMat :: [Sigma] -> [[Bool]]
    acommMat gs = [[acommQ g g' | g <- gs] | g' <- gs] -- TODO use asym for speedup
    -- is g cut by the region [i..i+l-1] ?
    cutByRegion :: Int -> Int -> Sigma -> Bool
    cutByRegion l i g | i+l<=n    = maybe False ((<i+l) . fst) (IntMap.lookupGE i g)
                                    && ( (fst $ IntMap.findMin g) < i || i+l <= (fst $ IntMap.findMax g) )
                      | otherwise = cutByRegion (n-l) (i+l-n) g
    -- local stabilizers cut by [i..i+n/2-1] where 0 <= i < n/2 due to symmetry
    -- TODO: I think this works even if n is odd, but I'm not certain
    cutStabs' :: IntMap [Sigma]
    cutStabs' = IntMap.unionsWith (flip (++)) $ map (\g -> IntMap.fromListWith (error "cutStabs'") $ map (,[g]) $ cutRegions $ IntMap.keys g) $ localStabs
      where cutRegions :: [Int] -> [Int]
            cutRegions is = concat $ map region0 $ zip is (tail is ++ [head is])
            region0 (i,j) | mod' (j-i) <= n_2 = (region `on` mod') (i+1)     (j+1)
                          | otherwise         = (region `on` mod') (j-n_2+1) (i+n_2+1)
            n_2 = div n 2
    region :: Int -> Int -> [Int]
    region i j = i<=j ? [i'..j'-1] $ [0..j'-1] ++ [i'..n-1]
      where i' = -lY * div (-i) lY
            j' =  lY * div   j  lY
    -- local and nonlocal stabilizers
    localStabs, nonlocalStabs :: [Sigma]
    (localStabs, nonlocalStabs) = partition localQ $ Map.keys $ gc'Ham stab
    localQ = local'Sigma n -- TODO get rid of this by using nlgc'Ham above?

-- eeAvg :: [[(Int,Double,Double)]] -> [(Int,Double,Double)]
-- eeAvg ees = [uncurry (fst3 $ head ees0,,) $ meanError $ map snd3 ees0 | ees0 <- transpose ees] -- |

randomBetas :: Double -> Int -> [F]
randomBetas gamma seed = randomBetas0 s0
  where s0 = ST.runST $ do
          gen <- Random.initialize $ Vector.fromList [fromIntegral seed]
          Random.save gen
        randomBetas0 s = r' : randomBetas0 s''
          where (r',s'') = ST.runST $ do 
                  gen <- Random.restore s
                  r   <- Random.beta alpha 1 gen
                  s'  <- Random.save gen
                  return (toF r, s')
                alpha = recip gamma

i_xs :: [Int] -> [Int] -> Int
i_xs ls xs = foldl' (\i (l,x) -> i*l + mod x l) 0 $ zip ls xs

xs_i :: [Int] -> Int -> [Int]
xs_i ls i0 = map snd $ init $ scanr (\l (i,_) -> divMod i l) (i0,0) ls

-- site -> randoms -> (SigmaTerms, unused_randoms)
type ModelGen = [Int] -> [F] -> ([([([Int],Int)],F)], [F])

init_generic'RG :: Double -> Int -> ([Int],ModelGen) -> RG
init_generic'RG gamma seed (ls,model) = init'RG $ fromLists'Ham ls $ mapMaybe f $ concat
                                      $ flip State.evalState (randomBetas gamma seed) $ mapM (State.state . model) $ mapM (\l -> [0..l-1]) ls
  where f :: ([([Int],Int)],F) -> Maybe ([(Int,Int)],F)
        f (_,0) = Nothing
        f (g,c) = Just (flip map g $ mapFst $ i_xs ls, c)

data Model = RandFerm | Ising | XYZ | MajChain | ToricCode -- |
  deriving (Eq, Show, Read)

model_gen :: Model -> [Int] -> [F] -> ([Int],ModelGen)
model_gen RandFerm ls [p] = (ls, gen)
  where n          = product ls
        gen [0] rs = foldl' gen' ([],rs) $
                       [ [([x],3)] | x <- [0..n-1] ] ++ -- |
                       [ [([x1],k1),([x2],k2)] ++ [([x],3) | x <- [x1+1..x2-1]] |
                         x1 <- [0..n-1], x2 <- [x1+1..n-1], k1 <- [1,2], k2 <- [1,2]]
        gen [_] rs = ([], rs)
        gen  _  _  = error "model_gen RandFerm"
        gen' (terms,rp:rh:rs') g = (if' (rp<p) [(g,rh)] [] ++ terms, rs') -- TODO p is affected by gamma
        gen' _                 _ = undefined
model_gen Ising ls [j,k,h] = (ls, gen)
  where (kj,kh) = (1,3)
        gen [x] (rj:rk:rh:rs) =
          ([ ([([x],kj),([x+1],kj)],j*rj), ([([x],kh),([x+1],kh)],k*rk), ([([x],kh)],h*rh) ], rs)
        gen [x,y] (rjx:rjy:rkx:rky:rh:rs) =
          ([ ([([x,y],kj),([x+1,y  ],kj)],j*rjx), ([([x,y],kh),([x+1,y  ],kh)],k*rkx), ([([x,y],kh)],h*rh),
             ([([x,y],kj),([x  ,y+1],kj)],j*rjy), ([([x,y],kh),([x  ,y+1],kh)],k*rky) ], rs)
        gen _ _ = error "model_gen Ising"
model_gen XYZ ls [jx,jy,jz] = (ls, gen)
  where gen [x] (rx:ry:rz:rs) =
          ([ ([([x],1),([x+1],1)],jx*rx), ([([x],2),([x+1],2)],jy*ry), ([([x],3),([x+1],3)],jz*rz) ], rs)
        gen _ _ = error "model_gen XYZ"
model_gen MajChain [lx] couplings = model_gen MajChain [lx,1] couplings
model_gen MajChain ls [t,t',h] = (ls++[3], gen)
  where (k1,k3, kh) = (1,3, 3)
        gen [i,a,0] (rt1:rt2:rt'1:rt'2:rh1:rh2:rs) =
            ([ (             ik1,t *rt1 ), (             ik2, t *rt2 ),   ([([i,a,1],kh)],-h*rh1),
               (([i,a,1],kh):ik1,t'*rt'1), (([i,a,2],kh):ik2,-t'*rt'2),   ([([i,a,2],kh)],-h*rh2) ], rs)
          where ik1 = [([i,a,0],k3)]
                ik2 = [([i,a,0],k1),([i+1,a,0],k1)]
        gen [_,7,_] _  = error "model_gen MajChain 7"
        gen [_,_,_] rs = ([], rs)
        gen _ _ = error "model_gen MajChain"
model_gen ToricCode ls [a,b,a',b'] = (ls++[2], gen)
  where (ka,kb) = (1,2)
        gen [x,y,0] (ra:rb:ra'1:ra'2:rb'1:rb'2:rs) =
          ([ ([([x,y,1],ka),([x+1,y,2],ka),([x,y+1,1],ka),([x,y,2],ka)],a*ra), ([([x,y,1],ka)],a'*ra'1), ([([x,y,2],ka)],a'*ra'2),
             ([([x,y,1],kb),([x-1,y,1],kb),([x,y-1,2],kb),([x,y,2],kb)],b*rb), ([([x,y,1],kb)],b'*rb'1), ([([x,y,2],kb)],b'*rb'2)], rs)
        gen [_,_,1] rs = ([],rs)
        gen _ _ = error "model_gen ToricCode"
model_gen _ _ _ = error "model_gen"

main :: IO ()
main = do
  (seed,model,ln2_ls,couplings,gamma) <- getArgs5 :: IO (Int, Model, [Int], [F], Double)
  let ls0      = map (2^) ln2_ls
      n        = n'RG rg0
      n_2      = div n 2
      rg0      = init_generic'RG gamma seed $ model_gen model ls0 couplings
      (rgs,rg) = mapSnd (snd . head) $ span ((>0) . fst) $ flip iterate (n,rg0)
               $ \(_n,_rg) -> let _n' = div _n 2 in (_n', nest (_n-_n') rgStep _rg)
    --rg       = runRG rg0
      xs       = map (2^) [0..head ln2_ls-1]
  
  putStr "model:     "; print $ show model
  putStr "Ls:        "; print ls0; 0 < n ? return() $ error "ln2_ls"
  putStr "couplings: "; print couplings
  putStr "Γ:         "; print gamma
  putStr "seed:      "; print seed
  putStr "max terms: "; print $ max_terms'RG rg
  
  putStr "entanglement entropy: " -- [(region size, entanglement entropy, error)]
  print $ ee1d xs $ stab'RG rg
  
  putStr "anderson correlator: " -- [(distance, [xyz -> (correlator,error)])]
  print $ [(x, [anderson_corr x (IntMap.singleton 0 k) $ stab'RG rg | k <- [1..3]]) | x <- xs]
  
  putStr "mid-RG coefficients: " -- [(remaining dimension,[coefficient])]
  print $ [(_n, on (++) Map.elems (gc'Ham $ ham'RG _rg) (diag'RG _rg)) | (_n,_rg) <- take 2 $ rgs++[(0,rg)]] -- |
  
  putStr "stabilizer length: " -- [(stabilizer length, coefficient)]
  print $ map (mapFst $ (\is -> (n-) $ foldl1' max $ map (\(i,j) -> abs $ mod (j-i-1) n) $ zip is $ tail is ++ [head is]) . IntMap.keys)
        $ Map.toList $ gc'Ham $ stab'RG rg
  
--putStr "stabilizer size: " -- [(stabilizer size, coefficient)]
--print $ map (mapFst IntMap.size) $ Map.toList $ gc'Ham $ stab'RG rg
  
  putStr "trash size: " -- [(typical sigma size, typical coefficient)]]
  print $ flip map (trash'RG rg) $ mapBoth (sqrt . sum . map (^(2::Int))) . mapFst (map $ fromIntegral . IntMap.size) . unzip
  
  let (d_diag, d_offdiag) = mapBoth (map $ mapFst IntMap.size) $ partition (all ((==3) . snd) . IntMap.toList . fst)
                          $ Map.toList $ gc'Ham $ c4s'Ham 1 (reverse $ g4s'RG $ rg) $ ham'RG rg0
  
  putStr "diagonalized offdiag: " -- [(size, coefficient)]
  print $ d_offdiag
  
  putStr "diagonalized diag: " -- [(size, coefficient)]
  print $ d_diag
  
--putStr "_unusedIs'RG: "
--print $ _unusedIs'RG rg
  
  putStr "small stabilizers: " 
  print $ take 20 $ map snd $ stab0'RG rg
  
--putStr "stabilizers: " 
--print $ stab'RG rg
  
  putStr "CPU time: "
  cpu_time <- getCPUTime
  print $ (1e-12::Double) * fromIntegral cpu_time

-- :main "0" "Ising" "[2]" "[2,1,1]" "1"
-- grep -v -e 'stabilizer length' -e 'anderson correlator' -e 'trash size'

-- troublesome:
-- ./SBRG 112 Ising [6] [2,1,1] 1
-- ./SBRG   3 Ising [9] [2,1,1] 1
