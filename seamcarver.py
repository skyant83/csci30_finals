#!/usr/bin/env python3

from picture import Picture
import math
import sys
from PIL import Image


class SeamCarver(Picture):
    ## TO-DO: fill in the methods below
    def energy(self, i: int, j: int) -> float:
        '''
        Return the energy of pixel at column i and row j
        '''
        if not (0 <= i < self.width() and 0 <= j < self.height()):
            raise IndexError("Pixel indices out of range.")

        left = self[i-1 if i > 0 else self.width() - 1, j]
        right = self[i+1 if i < self.width() - 1 else 0, j]
        up = self[i, j-1 if j > 0 else self.height() - 1] 
        down = self[i, j+1 if j < self.height() - 1 else 0] 

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
        dp = [[sys.maxsize] * width for _ in range(height)]
        path = [[0] * width for _ in range(height)]

        # Initialize the top row
        for i in range(width):
            dp[0][i] = energy_map[0][i]

        # Populate the DP table
        for j in range(1, height):
            for i in range(width):
                min_energy = dp[j - 1][i]
                if i > 0:
                    min_energy = min(min_energy, dp[j - 1][i - 1])
                if i < width - 1:
                    min_energy = min(min_energy, dp[j - 1][i + 1])

                dp[j][i] = energy_map[j][i] + min_energy

                # Track the path
                if min_energy == dp[j - 1][i]:
                    path[j][i] = i
                elif i > 0 and min_energy == dp[j - 1][i - 1]:
                    path[j][i] = i - 1
                else:
                    path[j][i] = i + 1

        # Backtrack to find the seam
        min_index = dp[-1].index(min(dp[-1]))
        seam = [0] * height
        for j in range(height - 1, -1, -1):
            seam[j] = min_index
            min_index = path[j][min_index]

        return seam

    def find_horizontal_seam(self) -> list[int]:
        '''
        Return a sequence of indices representing the lowest-energy horizontal seam
        '''
        # Transpose the image
        transposed_image = self.picture().transpose(Image.Transpose.TRANSPOSE)
        # Create a SeamCarver instance for the transposed image
        transposed_seam_carver = SeamCarver(transposed_image)
        # Find the vertical seam in the transposed image
        transposed_seam = transposed_seam_carver.find_vertical_seam()
        # Map the transposed seam back to horizontal seam
        horizontal_seam = transposed_seam
        
        return horizontal_seam

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
        '''
        Remove a horizontal seam from the picture
        '''
        # Transpose the image to swap width and height
        transposed_image = self.picture().transpose(Image.Transpose.TRANSPOSE)
        # Create a SeamCarver instance for the transposed image
        transposed_seam_carver = SeamCarver(transposed_image)
        # Remove the seam by utilizing the remove_vertical_seam
        transposed_seam_carver.remove_vertical_seam(seam)
        # Transpose the image back to its original dimensions
        updated_image = transposed_seam_carver.picture().transpose(Image.Transpose.TRANSPOSE)

        # Update the current object's pixel data directly
        self._width, self._height = updated_image.size
        self.clear()  # Clear the current pixel data
        for j in range(self._height):
            for i in range(self._width):
                self[i, j] = updated_image.getpixel((i, j))

class SeamError(Exception):
    pass
