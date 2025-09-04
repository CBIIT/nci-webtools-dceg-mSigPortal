#!/bin/bash

# Get git tag version
if git describe --tags --abbrev=0 >/dev/null 2>&1; then
    LATEST_TAG=$(git describe --tags --abbrev=0)
    VERSION=$(echo "$LATEST_TAG" | sed 's/msigportal_//g' | sed 's/_[0-9]\{8\}$//g')
    
    # Get the date of the latest tag
    TAG_DATE=$(git show --format="%cd" --date=format:"%B %d, %Y" --no-patch "$LATEST_TAG" 2>/dev/null || echo "Unknown")
else
    # Fallback if no tags
    VERSION="1.1.0-dev"
    TAG_DATE=$(date +"%B %d, %Y")
fi

# Create .env file with version information
cat > .env << EOF
APP_PATH=/mutational-signatures

# Application version and deployment information (auto-generated)
VITE_APP_VERSION=$VERSION
VITE_APP_LAST_UPDATE=$TAG_DATE
EOF

echo "Updated .env with:"
echo "VITE_APP_VERSION=$VERSION"
echo "VITE_APP_LAST_UPDATE=$TAG_DATE"

# Run the build
npm run build
