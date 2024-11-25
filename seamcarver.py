#!/usr/bin/env python3

from picture import Picture
import math
import sys

class SeamCarver(Picture):
    ## TO-DO: fill in the methods below
    def energy(self, i: int, j: int) -> float:
        '''
        Return the energy of pixel at column i and row j
        '''
        if not (0 <= i < self.width() and 0 <= j < self.height()):
            raise IndexError("Pixel indices out of range.")

        left = self[(i - 1) % self.width() if i > 0 else self.width() - 2, j]
        right = self[i % self.width() if i < self.width() else 0, j]
        up = self[i, (j - 1) % self.height() if j > 0 else self.height() - 2]
        down = self[i, j % self.height() if j < self.height() else 0]

        # Compute gradients for RGB components
        dx_r, dx_g, dx_b = (right[0] - left[0], right[1] - left[1], right[2] - left[2])
        dy_r, dy_g, dy_b = (down[0] - up[0], down[1] - up[1], down[2] - up[2])

        delta_x2 = dx_r**2 + dx_g**2 + dx_b**2
        delta_y2 = dy_r**2 + dy_g**2 + dy_b**2

        return math.sqrt(delta_x2 + delta_y2)
    
    
    def find_vertical_seam(self) -> list[int]:
        '''
        Return a sequence of indices representing the lowest-energy vertical seam
        '''
        width, height = self.width(), self.height()
        energy_map = [[self.energy(i, j) for i in range(width)] for j in range(height)]
        dp = [[0] * width for _ in range(height)]

        # Populate the DP table
        for j in range(height):
            for i in range(width):
                if j == 0:
                    dp[j][i] = energy_map[j][i]
                else:
                    # Calculate minimum energy from neighbors
                    min_energy = dp[j - 1][i]
                    if i > 0:
                        min_energy = min(min_energy, dp[j - 1][i - 1])
                    if i < width - 1:
                        min_energy = min(min_energy, dp[j - 1][i + 1])

                    dp[j][i] = energy_map[j][i] + min_energy

        # Backtrack to find the seam
        seam = [0] * height
        seam[-1] = dp[-1].index(min(dp[-1]))  # Start at the bottom row

        for j in range(height - 2, -1, -1):
            i = seam[j + 1]
            # Find the parent pixel in the previous row with the minimum energy
            min_index = i
            if i > 0 and dp[j][i - 1] < dp[j][min_index]:
                min_index = i - 1
            if i < width - 1 and dp[j][i + 1] < dp[j][min_index]:
                min_index = i + 1
            seam[j] = min_index

        return seam

    def find_horizontal_seam(self) -> list[int]:
        '''
        Return a sequence of indices representing the lowest-energy horizontal seam.
        '''
        width, height = self.width(), self.height()
        energy_map = [[self.energy(i, j) for i in range(width)] for j in range(height)]
        dp = [[sys.maxsize] * height for _ in range(width)]
        path = [[0] * height for _ in range(width)]

        for j in range(height):
            dp[0][j] = energy_map[j][0]

        for i in range(1, width):
            for j in range(height):
                min_energy = dp[i - 1][j]
                if j > 0:
                    min_energy = min(min_energy, dp[i - 1][j - 1])
                if j < height - 1:
                    min_energy = min(min_energy, dp[i - 1][j + 1])

                dp[i][j] = energy_map[j][i] + min_energy

                if min_energy == dp[i - 1][j]:
                    path[i][j] = j
                elif j > 0 and min_energy == dp[i - 1][j - 1]:
                    path[i][j] = j - 1
                else:
                    path[i][j] = j + 1

        min_index = dp[-1].index(min(dp[-1]))
        seam = [0] * width
        for i in range(width - 1, -1, -1):
            seam[i] = min_index
            min_index = path[i][min_index]

        return seam


    def remove_vertical_seam(self, seam: list[int]):
        '''
        Remove a vertical seam from the picture
        '''
        if len(seam) != self.height():
            raise SeamError("Invalid seam length.")
        if self.width() <= 1:
            raise SeamError("Cannot remove seam from image with width <= 1.")

        for j in range(self.height()):
            if j > 0 and abs(seam[j] - seam[j - 1]) > 1:
                raise SeamError("Invalid seam: adjacent entries differ by more than 1.")
            for i in range(seam[j], self.width() - 1):
                self[i, j] = self[i + 1, j]
            del self[self.width() - 1, j]

        self._width -= 1


    def remove_horizontal_seam(self, seam: list[int]):
        if len(seam) != self.width():
            raise SeamError("Invalid seam length.")
        if self.height() <= 1:
            raise SeamError("Cannot remove seam from image with height <= 1.")

        for i in range(self.width()):
            if i > 0 and abs(seam[i] - seam[i - 1]) > 1:
                raise SeamError("Invalid seam: adjacent entries differ by more than 1.")
            for j in range(seam[i], self.height() - 1):
                self[i, j] = self[i, j + 1]
            del self[i, self.height() - 1]

        self._height -= 1


class SeamError(Exception):
    pass
