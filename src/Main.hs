{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE NoImplicitPrelude #-}
module Main where

import qualified Prelude as P
import Data.Array.Accelerate.Data.Complex as A
import Data.Array.Accelerate as A
-- import qualified Data.Array.Accelerate.Fourier.Adhoc as A
import qualified Data.Array.Accelerate.Interpreter as I
import Data.Array.Accelerate.Interpreter as I
import Data.Array.Accelerate.Data.Fold
import Data.Array.Accelerate.Data.Monoid
import qualified Data.List as L
import qualified Data.Array.Accelerate.LLVM.PTX as PTX

import qualified Data.Csv as CSV
import qualified Data.ByteString.Lazy as BS
import qualified Statistics.Sample as S
import qualified Data.Vector as V

type Arr a = Acc (Vector a)

mean :: (Floating a, FromIntegral Int a) => Arr a -> Acc (Scalar a)
mean xs = unit $ the (fold (+) 0 xs) / (fromIntegral $ length xs)

smooth :: (Floating a, FromIntegral Int a) => Arr a -> Exp Int -> Arr a
smooth xs res = error "LOL"

sumF :: (Elt a, Lift Exp a, P.Num (Exp a)) => Fold (Exp a) (Exp a)
sumF = Fold (lift . Sum) (getSum . unlift)

lengthF :: (Elt a, Lift Exp a, P.Num (Exp a)) => Fold (Exp a) (Exp a)
lengthF = Fold (\_ -> 1) (getSum . unlift)

meanF :: (Elt a, Lift Exp a, P.Fractional (Exp a)) => Fold (Exp a) (Exp a)
meanF = sumF / lengthF

valF :: (Elt a, Monoid (Exp a)) => Fold (Exp a) (Exp a)
valF = Fold P.id P.id

-- | E(X^2)
sumOfSquaresF :: (Elt a, Lift Exp a, P.Num (Exp a)) => Fold (Exp a) (Exp a)
sumOfSquaresF = Fold (lift . Sum . (^(2 :: Exp Int))) (getSum . unlift)

-- | sqrt (E(X^2) - E(X) ^ 2)
standardDeviationF :: (Elt a, Lift Exp a, P.Floating (Exp a)) => Fold (Exp a) (Exp a)
standardDeviationF = sqrt $ (sumOfSquaresF / lengthF) - (meanF * meanF)

-- | E(X - mu) ^ 4
-- differenceQuadsF :: (Elt a, Lift Exp a, P.Fractional (Exp a), Monoid (Exp a)) => Fold (Exp a) (Exp a)
-- differenceQuadsF = (\m -> Fold (\x -> lift $ Sum ((x - m) ^ e)) (getSum . unlift)) P.<$> meanF
--   where
--     e = 4 :: Exp Int

-- | E((X - u)^2)^2
-- otherThing :: (Elt a, Lift Exp a, P.Fractional (Exp a), Monoid (Exp a)) => Fold (Exp a) (Exp a)
-- otherThing = Fold (lift . Sum . (^(2 :: Exp Int))) (getSum . unlift)

-- kurtosisF :: (Elt a, Lift Exp a, P.Fractional (Exp a), Monoid (Exp a)) => Fold (Exp a) (Exp a)
-- kurtosisF = (valF - meanF)

-- | Kurt(X) = E(X - u)^4/(E((X - u)^2)^2)
-- (mean (map (\x -> x - u) xs)) ^ four / ((mean (map (\x -> (x - u) ^ two) xs))^two)
kurtosis :: (Floating a, Elt a, FromIntegral Int a) => Arr a -> Exp a
kurtosis xs = (u4 / (variance P.^ 2))
  where
    u4 = the $ compute $ mean $ map (\x -> (x - u) ^ (4 :: Exp Int)) xs
    variance = the $ compute $ mean $ map (\x -> (x - u) ^ (2 :: Exp Int)) xs
    u = the $ compute $ mean xs

diffs :: (Elt a, P.Num (Exp a)) => Arr a -> Arr a
diffs xs = zipWith (-) (tail xs) xs

roughness :: (Elt a, P.Floating (Exp a), FromIntegral Int a, Lift Exp a) => Arr a -> Acc (Scalar a)
roughness xs = compute $ runFold standardDeviationF $ diffs xs

-- smoothSimple :: (Elt a, Floating a, FromIntegral Int a, Ord a) => Arr a -> Exp Int -> Arr a
-- smoothSimple xs resolution =
--   ifThenElse (resolution < length xs)
--     smaLargeResolution
--     (sma bestWindowSize 1 xs)
--   where
--     -- THIS IS WRONG. IT NEEDS TO MINIMIZE THE ROUGHNESS
--     bestWindowSize =  the $ minimum $ map fst $ afst $ filter acceptableWindow $ map (go xs) ws
--     acceptableWindow b =
--       (r < originalRoughness) && (k >= originalKurt)
--       where
--         r = fst $ snd b
--         k = snd $ snd b
--     smaLargeResolution = sma (length xs `div` resolution) (length xs `div` resolution) xs
--     originalKurt = kurtosis xs
--     originalRoughness = the $ roughness xs
--     end :: Exp Int
--     end = length xs `div` 10
--     ws :: Arr Int
--     ws = enumFromN (index1 $ end - 2) 2

-- smoothSimple' :: (Elt a, Floating a, FromIntegral Int a, Ord a) => Arr a -> Exp Int -> Arr a
-- smoothSimple' xs resolution = sma bestWindowSize 1 xs
--   where
--     bestWindowSize =
--       fst $ while (\a -> fst a >= 2) godude $ lift ((length xs `div` resolution), (originalRoughness, 1 :: Exp Int))
--     originalKurt = kurtosis xs
--     originalRoughness = the $ roughness xs
--     godude a =
--       let
--         w = fst a
--         minObj = fst $ snd a
--         smoothed = sma w 1 xs
--         r = the $ roughness smoothed
--         kurt = kurtosis smoothed
--       in ifThenElse (r < minObj && kurt >= originalKurt) (lift (w-1, (r, w))) (lift (w-1, snd a))

-- go :: (Elt a, P.Floating (Exp a), FromIntegral Int a) => Arr a -> Exp Int -> Exp (Int, (a, a))
-- go xs w = lift (w,(the $ roughness smoothed, kurtosis smoothed))
--   where
--     smoothed = sma w 1 xs

sma :: (Floating a, FromIntegral Int a) => Exp Int -> Exp Int -> Arr a -> Arr a
sma range slide xs = mavg range xs
  where
    skip ys =
      map snd. afst
        $ filter (\y -> fst y `mod` slide == 0)
          $ imap (\i y -> lift (unindex1 i, y)) ys

computeStuff :: (Floating a, FromIntegral Int a, Lift Exp a) => Exp Int -> Exp Int -> Arr a -> Exp (a, a)
computeStuff resolution window xs = lift (the $ roughness smoothed, kurtosis smoothed)
  where
    smoothed = sma window resolution xs

kurtAndRoughness :: (Floating a, FromIntegral Int a, Lift Exp a) => Arr a -> Exp (a, a)
kurtAndRoughness xs = lift (the $ roughness xs, kurtosis xs)

mavg :: (Floating a, FromIntegral Int a) => Exp Int -> Arr a -> Arr a
mavg k xs =
  map (/ fromIntegral k)
    $ scanl (+) (the $ sum h)
      $ zipWith (-) t xs
  where
    h = take k xs
    t = drop k xs

-- transform :: Arr (Complex Double) -> Arr (Complex Double)
-- transform = A.transform A.forward

-- inverseTransform :: Arr (Complex Double) -> Arr (Complex Double)
-- inverseTransform = A.transform A.inverse

binarySearch :: Int -> Int -> Arr (Complex Double) -> a -> a -> a -> ()
binarySearch h t xs minObj originalKurt windowSize = ()

main :: P.IO ()
main = do
  P.Right xs <- CSV.decode CSV.HasHeader P.<$> BS.readFile "data.csv" :: P.IO (P.Either P.String (V.Vector (P.String, Double)))

  let
    -- use accelerate-io
    v = V.toList $ P.snd P.<$> xs :: [Double]
    vec = use $ fromList (Z :. (V.length xs)) v
    ws = P.reverse [2..V.length xs `div` 10]
    [(minObj, originalKurt)] = toList $ PTX.runN $ unit $ kurtAndRoughness vec
    runOnWindow w = P.head $ toList $ PTX.runN $ unit $ computeStuff 1 (constant w) vec
    roughsAndKurts = (\w -> (w, runOnWindow w)) P.<$> ws
    bestWindow@(_, w) = L.foldl' goman (minObj, 1) roughsAndKurts
    goman old@(obj, _) (w,(r, k)) = if (r P.< obj P.&& k P.>= originalKurt) then (r,w) else old
    smoothed = PTX.run $ sma (constant w) 1 vec
  P.print $ L.take 100 roughsAndKurts
  P.print $ (minObj, originalKurt)
  P.print $ V.length xs
  P.print smoothed
  -- P.print $ S.kurtosis $ (V.fromList v :: V.Vector Double)
    -- bestWindowSize =  the $ minimum $ map fst $ afst $ filter acceptableWindow $ map (go xs) ws
    -- acceptableWindow b =
    --   (r < originalRoughness) && (k >= originalKurt)
    --   where
    --     r = fst $ snd b
    --     k = snd $ snd b
    -- smaLargeResolution = sma (length xs `div` resolution) (length xs `div` resolution) xs
    -- originalKurt = kurtosis xs
    -- originalRoughness = the $ roughness xs
    -- end :: Exp Int
    -- end = length xs `div` 10
    -- ws :: Arr Int
    -- ws = enumFromN (index1 $ end - 2) 2

mavg' :: P.Fractional b => Int -> [b] -> [b]
mavg' k lst = P.map (/ P.fromIntegral k) $ L.scanl' (+) (P.sum h) $ P.zipWith (-) t lst
  where (h, t) = L.splitAt k lst