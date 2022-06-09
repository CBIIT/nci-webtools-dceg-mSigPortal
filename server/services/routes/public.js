const { Router } = require('express');
const cors = require('cors');

const router = Router();

router.use(cors());
router.get('/ping', (req, res) => res.send(true));

module.exports = router;
